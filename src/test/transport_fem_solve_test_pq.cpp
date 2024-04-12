#include "cfd24_test.hpp"
#include "cfd24/grid/unstructured_grid2d.hpp"
#include "cfd24/grid/grid1d.hpp"
#include "cfd24/mat/csrmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/mat/lodmat.hpp"
#include "utils/filesystem.hpp"
#include "cfd24/fem/fem_assembler.hpp"
#include "cfd24/fem/elem1d/segment_linear.hpp"
#include "cfd24/fem/elem2d/triangle_linear.hpp"
#include "cfd24/fem/elem2d/quadrangle_linear.hpp"
#include "cfd24/numeric_integration/segment_quadrature.hpp"
#include "cfd24/numeric_integration/square_quadrature.hpp"
#include "cfd24/numeric_integration/triangle_quadrature.hpp"
#include "cfd24/fem/fem_numeric_integrals.hpp"
#include "cfd24/debug/printer.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include <algorithm>
static double limiter(double r) {
	int number_limited = 3;
	switch (number_limited)
	{
		// upwind
		case 1: return 0.0;
		// sym
		case 2: return 1.0;
		// minmod
		case 3: return std::max(0.0, std::min(r, 1.0));
		// superbee
		case 4: return std::max(0.0, std::max(std::min(2.0, std::abs(r)), std::min(1.0, 2.0 * std::abs(r))));
		default: throw std::runtime_error("Invalid limiter");
	}
}
using namespace cfd;

class ATestTransport2FemWorker {
public:
	struct Edge {
		size_t i, j;
		size_t ii_addr, jj_addr, ij_addr, ji_addr;
		void reverse() {
			std::swap(i, j);
			std::swap(ii_addr, jj_addr);
			std::swap(ij_addr, ji_addr);
		}
	};

	virtual double init_solution(Point p) const {
		double x = p.x();
		double y = p.y();
		return (x >= 0.05 && x <= 0.15 && y >= 0.05 && y <= 0.15) ? 1.0 : 0.0;
		//return (x >= 0.1 && x <= 0.2) ? 1.0 : 0.0;
	}
	virtual double exact_solution(Point p) const {
		return init_solution(p - _time * Point(1, 1));
		//return init_solution(p - _time * Point(1, 0));
	}

	ATestTransport2FemWorker(const IGrid& g) : _grid(g), _fem(build_fem(g)), _u(_grid.n_points()){};

	void initialize(){
		// 1. Lumped mass
		CsrMatrix full_mass(_fem.stencil());
		for (size_t ielem = 0; ielem < _fem.n_elements(); ++ielem) {
			const FemElement& elem = _fem.element(ielem);
			std::vector<double> local_mass = elem.integrals->mass_matrix();
			_fem.add_to_global_matrix(ielem, local_mass, full_mass.vals());
		}
		_mass = full_mass.mult_vec(std::vector<double>(_grid.n_points(), 1));
		// 2. 2nd order transport
		std::vector<double> vx(_grid.n_points(), 1.0);
		std::vector<double> vy(_grid.n_points(), 1.0);
		_K = _fem.stencil();
		for (size_t ielem = 0; ielem < _fem.n_elements(); ++ielem) {
			const FemElement& elem = _fem.element(ielem);
			std::vector<double> local_transport = elem.integrals->transport_matrix(
				_fem.local_vector(ielem, vx), _fem.local_vector(ielem, vy));
			_fem.add_to_global_matrix(ielem, local_transport, _K.vals());
		}
		for (double& v : _K.vals()) v *= -1;
		// 3. edges
		const auto& addr = _K.addr();
		const auto& cols = _K.cols();
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			for (size_t a = addr[i]; a < addr[i + 1]; ++a) {
				size_t j = cols[a];
				if (j >= i) break;
				_edges.push_back(Edge{
					i, j,
					_K.get_address(i, i), _K.get_address(j, j),
					a, _K.get_address(j, i)
					});
			}
		}
		// 4. upwind transport + edge direction
		_L = _K;
		for (auto& edge : _edges) {
			double& lii = _L.vals()[edge.ii_addr];
			double& ljj = _L.vals()[edge.jj_addr];
			double& lij = _L.vals()[edge.ij_addr];
			double& lji = _L.vals()[edge.ji_addr];
			double dij = std::max(0.0, std::max(-lij, -lji));
			lii -= dij;
			lij += dij;
			ljj -= dij;
			lji += dij;
			if (lji < lij) edge.reverse();
		}

		// 5. initial solution
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			_u[i] = init_solution(_grid.point(i));
		}
	}

	void step(double tau) {
		_time += tau;
		impl_step(tau);
	}

	void save_vtk(const std::string& filename) {
		// save grid
		_grid.save_vtk(filename);

		// save numerical solution
		VtkUtils::add_point_data(_u, "numerical", filename);

		// save exact solution
		std::vector<double> exact(_grid.n_points());
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			exact[i] = exact_solution(_grid.point(i));
		}
		VtkUtils::add_point_data(exact, "exact", filename);
	}

	double current_time() const {
		return _time;
	}

	double compute_norm2() const {
		// sum
		double sum = 0;
		double vol = 0;
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			double diff = _u[i] - exact_solution(_grid.point(i));
			sum += diff * diff * _mass[i];
			vol += _mass[i];
		}
		return std::sqrt(sum / vol);
	}

	double compute_norm_max() const {
		return 1.0 - *std::max_element(_u.begin(), _u.end());
	};

	double recommended_tau(double theta) const {
		std::vector<double> diag = _L.diagonal();
		double minv = 0;
		for (size_t i = 0; i < diag.size(); ++i) {
			minv = std::min(minv, diag[i] / _mass[i]);
		}
		return -1.0 / (1 - theta) / minv;
	};
protected:
	const IGrid& _grid;
	FemAssembler _fem;
	std::vector<double> _u;
	std::vector<double> _mass;
	CsrMatrix _K, _L, _K_star;
	std::vector<Edge> _edges;
	double _time = 0;

private:
	static FemAssembler build_fem(const IGrid& grid) {
		size_t n_bases = grid.n_points();
		std::vector<FemElement> elements;
		std::vector<std::vector<size_t>> tab_elem_basis;

		// elements
		for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
			std::vector<size_t> ipoints = grid.tab_cell_point(icell);
			tab_elem_basis.push_back(ipoints);
			if (ipoints.size() == 2) {
				Point p0 = grid.point(ipoints[0]);
				Point p1 = grid.point(ipoints[1]);

				auto geom = std::make_shared<SegmentLinearGeometry>(p0, p1);
				auto basis = std::make_shared<SegmentLinearBasis>();
				const Quadrature* quadrature = quadrature_segment_gauss2();
				auto integrals = std::make_shared<NumericElementIntegrals>(quadrature, geom, basis);
				elements.push_back(FemElement{ geom, basis, integrals });
			}
			else if (ipoints.size() == 3) {
				Point p0 = grid.point(ipoints[0]);
				Point p1 = grid.point(ipoints[1]);
				Point p2 = grid.point(ipoints[2]);

				auto geom = std::make_shared<TriangleLinearGeometry>(p0, p1, p2);
				auto basis = std::make_shared<TriangleLinearBasis>();
				const Quadrature* quadrature = quadrature_triangle_gauss2();
				auto integrals = std::make_shared<NumericElementIntegrals>(quadrature, geom, basis);
				elements.push_back(FemElement{ geom, basis, integrals });
			}
			else if (ipoints.size() == 4) {
				Point p0 = grid.point(ipoints[0]);
				Point p1 = grid.point(ipoints[1]);
				Point p2 = grid.point(ipoints[2]);
				Point p3 = grid.point(ipoints[3]);

				auto geom = std::make_shared<QuadrangleLinearGeometry>(p0, p1, p2, p3);
				auto basis = std::make_shared<QuadrangleLinearBasis>();
				const Quadrature* quadrature = quadrature_square_gauss2();
				auto integrals = std::make_shared<NumericElementIntegrals>(quadrature, geom, basis);
				elements.push_back(FemElement{ geom, basis, integrals });
			}
			else {
				throw std::runtime_error("invalid fem grid");
			}
		}

		return FemAssembler(n_bases, elements, tab_elem_basis);
	}
	virtual void impl_step(double tau) = 0;
};

///////////////////////////////////////////////////////////////////////////////
// Explicit transport solver
///////////////////////////////////////////////////////////////////////////////

class TestTransport2FemUpwindExplicit : public ATestTransport2FemWorker {
public:
	TestTransport2FemUpwindExplicit(const IGrid& g) : ATestTransport2FemWorker(g) {}
private:
	void impl_step(double tau) override {
		std::vector<double> uold = _u;
		std::vector<double> P_p(_grid.n_points(), 0.0);
		std::vector<double> P_m(_grid.n_points(), 0.0);
		std::vector<double> Q_p(_grid.n_points(), 0.0);
		std::vector<double> Q_m(_grid.n_points(), 0.0);
		for (auto& edge : _edges) {
			double du = _u[edge.j] - _u[edge.i];
			double kij = _K.vals()[edge.ij_addr];
			P_p[edge.i] += std::min(0.0, kij) * std::min(0.0, du);
			P_m[edge.i] += std::min(0.0, kij) * std::max(0.0, du);
			Q_m[edge.i] += std::max(0.0, kij) * std::min(0.0, du);
			Q_p[edge.i] += std::max(0.0, kij) * std::max(0.0, du);

			double kji = _K.vals()[edge.ji_addr];
			P_p[edge.j] += std::min(0.0, kji) * std::min(0.0, -du);
			P_m[edge.j] += std::min(0.0, kji) * std::max(0.0, -du);
			Q_m[edge.j] += std::max(0.0, kji) * std::min(0.0, -du);
			Q_p[edge.j] += std::max(0.0, kji) * std::max(0.0, -du);
		}
		CsrMatrix Fa(_fem.stencil());
		for (auto& edge: _edges) {
			double& lij = _K.vals()[edge.ij_addr];
			double& lji = _K.vals()[edge.ji_addr];
			double dij = std::max(0.0, std::max(-lij, -lji));
			if (lji < lij) edge.reverse();
			//==============================================================->
			double r_i;
			double fa_ji;
			double fa_ij;
			//auto x = P_p[edge.i];
			//auto y = P_m[edge.i];
			if (_u[edge.i] > _u[edge.j]){ 
				double denum = P_p[edge.i];
				if (denum == 0) denum = -1e-6;
				r_i = Q_p[edge.i] / denum;
			} else {
				double denum = P_m[edge.i];
				if (denum == 0) denum = -1e-6;
				r_i = Q_m[edge.i] / denum;
			}

			fa_ij = -std::min(limiter(r_i) * dij, lji);
			fa_ji = -fa_ij;
			Fa.vals()[edge.ij_addr] = fa_ij;
			Fa.vals()[edge.ii_addr] -= fa_ij;
			Fa.vals()[edge.ji_addr] = fa_ji;
			Fa.vals()[edge.jj_addr] -= fa_ji;
			//=================================================================
		}
		_K_star = _fem.stencil();
		for (size_t a = 0; a < _K_star.n_nonzeros(); ++a){
			//_K_star.vals()[a] = _K.vals()[a];
			_K_star.vals()[a] = _L.vals()[a] + Fa.vals()[a];

		}
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			_u[i] += tau / _mass[i] * _K_star.mult_vec(i, uold);
		}
	}
};

TEST_CASE("Transport 2D fem solver, explicit", "[transport2-fem-upwind-explicit]") {
	std::cout << std::endl << "--- cfd24_test [transport2-fem-upwind-explicit] --- " << std::endl;
	double tend = 0.5;

	// solver
	std::string fn = tmp_directory_file("tetragrid_10k.vtk");
	//std::string fn = test_directory_file("tetragrid_500.vtk");
	auto g = UnstructuredGrid2D::vtk_read(fn);
	TestTransport2FemUpwindExplicit worker(g);
	worker.initialize();
	double tau = worker.recommended_tau(0.0);

	// saver
	VtkUtils::TimeSeriesWriter writer("transport2-fem-upwind-explicit");
	std::string out_filename = writer.add(0);
	worker.save_vtk(out_filename);

	double n2 = 0, nm = 0;
	while (worker.current_time() < tend - 1e-6) {
		// solve problem
		worker.step(tau);
		// export solution to vtk
		out_filename = writer.add(worker.current_time());
		worker.save_vtk(out_filename);

		n2 = worker.compute_norm2();
		nm = worker.compute_norm_max();

		std::cout << worker.current_time() << " " << n2 << " " << nm << std::endl;
	};
	CHECK(nm == Approx(0.8961).margin(1e-4));
}


TEST_CASE("tvd-1d", "[tvd-ex]"){
	auto g = Grid1D(0, 1, 20);
	double h = g.cell_volume(0);
	double tau = h/3.0;

	auto exact = [](double x, double t)->double{
		return (x-t >= 0.1 && x-t <= 0.2) ? 1.0 : 0.0;
	};
	
	double tm = 0;
	std::vector<double> u(g.n_points(), 0);
	for (size_t i=0; i<g.n_points(); ++i){
		u[i] = exact(g.point(i).x(), 0);
	}
	std::vector<double> Fp(g.n_points(), 0);

	// saver
	VtkUtils::TimeSeriesWriter writer("tvd");
	auto saver = [&](double t){
		std::string out_filename = writer.add(t);
		g.save_vtk(out_filename);
		VtkUtils::add_point_data(u, "numerical", out_filename);
		std::vector<double> u_exact(g.n_points());
		for (size_t i = 0; i < g.n_points(); ++i) {
			u_exact[i] = exact(g.point(i).x(), t);
		}
		VtkUtils::add_point_data(u_exact, "exact", out_filename);
	};
	saver(0);

	while (tm < 0.5 - 1e-6){
		for (size_t i=0; i<g.n_points()-1; ++i){
			double denum = u[i+1] - u[i];
			if (denum == 0) denum = 1e-6;
			double r = (i == 0) ? 1 : (u[i] - u[i-1])/denum;
			Fp[i] = u[i] + 0.5*limiter(r)*(u[i+1] - u[i]);
		}
		for (size_t i=1; i<g.n_points()-1; ++i){
			u[i] -= tau*(Fp[i] - Fp[i-1])/h;
		}
		tm += tau;
		saver(tm);
		std::cout << tm << std::endl;
	}
}

class T2: public TestTransport2FemUpwindExplicit{
public:
	T2(const IGrid& g): TestTransport2FemUpwindExplicit(g){};

	double init_solution(Point p) const override {
		double x = p.x();
		return (x >= 0.1 && x <= 0.2) ? 1.0 : 0.0;
	}
	double exact_solution(Point p) const override {
		return init_solution(p - _time * Vector(1));
	}
};

TEST_CASE("tvd-1d-fem", "[tvd-ex-fem]"){
	double tend = 0.5;

	// solver
	auto g = Grid1D(0, 1, 20);
	T2 worker(g);
	worker.initialize();
	double tau = g.cell_volume(0)/3.0;

	// saver
	VtkUtils::TimeSeriesWriter writer("tvd2");
	std::string out_filename = writer.add(0);
	worker.save_vtk(out_filename);

	double n2 = 0, nm = 0;
	while (worker.current_time() < tend - 1e-6) {
		// solve problem
		worker.step(tau);
		// export solution to vtk
		out_filename = writer.add(worker.current_time());
		worker.save_vtk(out_filename);

		n2 = worker.compute_norm2();
		nm = worker.compute_norm_max();

		std::cout << worker.current_time() << " " << n2 << " " << nm << std::endl;
	};
}

TEST_CASE("", "[a]"){
	RegularGrid2D grid(0, 1, 0, 1, 100, 100);
	double eps = 1e-3;
	double t0 = 0.2;

	auto exact = [&](Point p, double t)->double{
		double x = p.x();
		double y = p.y();

		double x0 = 0.8;
		double y0 = 0.5;

		double x1 = 0.5 + (x0 - 0.5) * cos(t) - (y0 - 0.5) * sin(t);
		double y1 = 0.5 + (x0 - 0.5) * sin(t) + (y0 - 0.5) * cos(t);

		double r2 = (x-x1)*(x-x1) + (y-y1)*(y-y1);
		return exp(-r2/(4*eps*(t+t0)))/(4*M_PI*eps*(t+t0));
	};

	auto diff = [&](Point p, double t)->double{
		constexpr double tau = 1e-6;
		constexpr double h = 1e-6;

		double u = exact(p, t);
		double u_t = exact(p, t+tau);
		double u_xp = exact(p + Point(h, 0), t);
		double u_yp = exact(p + Point(0, h), t);
		double u_xm = exact(p - Point(h, 0), t);
		double u_ym = exact(p - Point(0, h), t);

		double vx = -p.y() + 0.5;
		double vy = +p.x() - 0.5;

		return (u_t - u)/tau
			+ vx * (u_xp - u_xm)/(2*h)
			+ vy * (u_yp - u_ym)/(2*h)
			+ eps/h/h*(4*u-u_xm-u_xp-u_ym-u_yp);
	};

	// saver
	VtkUtils::TimeSeriesWriter writer("tmp");
	for (double t = 0.0; t < 2*M_PI; t+=0.1){
		std::vector<double> u(grid.n_points());
		std::vector<double> d(grid.n_points());
		for (size_t i=0; i<grid.n_points(); ++i){
			u[i] = exact(grid.point(i), t);
			d[i] = diff(grid.point(i), t);
		}

		std::string filename = writer.add(t);
		grid.save_vtk(filename);
		VtkUtils::add_point_data(u, "exact", filename);
		VtkUtils::add_point_data(d, "diff", filename);
	}
}
