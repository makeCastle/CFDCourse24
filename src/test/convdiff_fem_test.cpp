#include "cfd24_test.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/grid/grid1d.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/mat/lodmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "utils/vecmat.hpp"
#include "cfd24/debug/tictoc.hpp"
#include "cfd24/debug/saver.hpp"
#include "cfd24/fvm/fvm_assembler.hpp"
#include "utils/filesystem.hpp"
#include "cfd24/grid/unstructured_grid2d.hpp"
#include "cfd24/fvm/fvm_dpdn_boundary.hpp"
#include "cfd24/fem/fem_assembler.hpp"
#include "cfd24/fem/fem_assembler.hpp"
#include "cfd24/fem/elem1d/segment_linear.hpp"
#include "cfd24/fem/elem1d/segment_quadratic.hpp"
#include "cfd24/fem/elem1d/segment_cubic.hpp"
#include "cfd24/fem/elem2d/triangle_linear.hpp"
#include "cfd24/fem/elem2d/triangle_quadratic.hpp"
#include "cfd24/fem/elem2d/quadrangle_linear.hpp"
#include "cfd24/fem/elem2d/quadrangle_quadratic.hpp"
#include "cfd24/numeric_integration/segment_quadrature.hpp"
#include "cfd24/numeric_integration/square_quadrature.hpp"
#include "cfd24/numeric_integration/triangle_quadrature.hpp"
#include "cfd24/fem/fem_numeric_integrals.hpp"
#include "cfd24/debug/printer.hpp"
#include "cfd24/fem/fem_sorted_cell_info.hpp"
#include<numeric>

using namespace cfd;

namespace {

FemAssembler build_fem1(unsigned power, const IGrid& grid){
	std::vector<FemElement> elements;
	std::vector<std::vector<size_t>> tab_elem_basis;

	// elements
	for (size_t icell=0; icell < grid.n_cells(); ++icell){
		std::vector<size_t> ipoints = grid.tab_cell_point(icell);
		if (ipoints.size() != 2){
			throw std::runtime_error("invalid fem grid");
		}
		tab_elem_basis.push_back(ipoints);

		Point p0 = grid.point(ipoints[0]);
		Point p1 = grid.point(ipoints[1]);
			
		std::shared_ptr<IElementBasis> basis;
		const Quadrature* quadrature;
		if (power == 1){
			basis = std::make_shared<SegmentLinearBasis>();
			quadrature = quadrature_segment_gauss2();
		} else if (power == 2){
			basis = std::make_shared<SegmentQuadraticBasis>();
			quadrature = quadrature_segment_gauss3();
			tab_elem_basis.back().push_back(grid.n_points() + icell);
		} else if (power == 3){
			basis = std::make_shared<SegmentCubicBasis>();
			quadrature = quadrature_segment_gauss4();
			tab_elem_basis.back().push_back(grid.n_points() + 2*icell);
			tab_elem_basis.back().push_back(grid.n_points() + 2*icell+1);
		}
		auto geom = std::make_shared<SegmentLinearGeometry>(p0, p1);
		auto integrals = std::make_shared<NumericElementIntegrals>(quadrature, geom, basis);
		elements.push_back(FemElement{geom, basis, integrals});
	}

	size_t n_bases = grid.n_points();
	if (power == 2){
		n_bases += grid.n_cells();
	} else if (power == 3){
		n_bases += 2*grid.n_cells();
	}
	return FemAssembler(n_bases, elements, tab_elem_basis);
}

}


TEST_CASE("Stationar 1d convection-diffustion", "[stat-convdiff-fem-supg]"){
	double norm2;
	double norm_max;
	//for (size_t n: {10, 20, 30, 40, 50, 70, 100, 200, 500, 1000})
	size_t n = 10;
	{
		double h = 1.0/n;
		Grid1D grid(0, 1, n);
		FemAssembler fa = build_fem1(1, grid);
		double eps = 1e-2;

		auto exact = [eps](Point p){
			double x = p.x();
			return x - std::expm1(x/eps)/std::expm1(1/eps);
		};

		std::vector<double> vx(fa.n_bases(), 1.0);
		CsrMatrix K = fa.stencil();
		CsrMatrix D = fa.stencil();
		CsrMatrix M = fa.stencil();

		CsrMatrix R = fa.stencil();
		std::vector<double> rf(fa.n_bases(), 0);

		for (size_t ielem=0; ielem < fa.n_elements(); ++ielem){
			const FemElement& elem = fa.element(ielem);
			auto v = fa.local_vector(ielem, vx);
			// transport
			std::vector<double> local_transport = elem.integrals->transport_matrix(v);
			fa.add_to_global_matrix(ielem, local_transport, K.vals());
			// diffusion
			std::vector<double> local_diff = elem.integrals->stiff_matrix();
			fa.add_to_global_matrix(ielem, local_diff, D.vals());
			// mass
			std::vector<double> local_mass = elem.integrals->mass_matrix();
			fa.add_to_global_matrix(ielem, local_mass, M.vals());
			// transport supg
			std::vector<double> local_supg = elem.integrals->transport_matrix_stab_supg(v);
			fa.add_to_global_matrix(ielem, local_supg, R.vals());
		}

		// lhs
		CsrMatrix A = fa.stencil();
		double s = 0.1*h;
		for (size_t a=0; a<A.n_nonzeros(); ++a){
			A.vals()[a] = K.vals()[a] + eps * D.vals()[a] + s*R.vals()[a];
		}

		// rhs
		std::vector<double> rhs = M.mult_vec(std::vector<double>(fa.n_bases(), 1.0));
		for (size_t i=0; i<rhs.size(); ++i){
			rhs[i] -= s*rf[i];
		}

		// Boundary condition
		A.set_unit_row(0);
		rhs[0] = exact(grid.point(0));
		A.set_unit_row(grid.n_points()-1);
		rhs[grid.n_points()-1] = exact(grid.point(grid.n_points()-1));

		// Solution
		std::vector<double> u;
		AmgcMatrixSolver::solve_slae(A, rhs, u);

		// exact
		std::vector<double> ex;
		for (size_t i=0; i<fa.n_bases(); ++i) ex.push_back(exact(fa.reference_point(i)));

		// saver
		std::string filename = "convdiff.vtk";
		grid.save_vtk(filename);
		VtkUtils::add_point_data(u, "numerical", filename, grid.n_points());
		VtkUtils::add_point_data(ex, "exact", filename, grid.n_points());

		std::vector<double> diff2(u.size());
		for (size_t i=0; i<diff2.size(); ++i) diff2[i] = (u[i] - ex[i])*(u[i] - ex[i]);
		diff2 = M.mult_vec(diff2);
		norm2 = std::sqrt(std::accumulate(diff2.begin(), diff2.end(), 0.0));
	
		norm_max = 0;
		for (size_t i=0; i<u.size(); ++i){
			norm_max = std::max(norm_max, std::abs(u[i] - ex[i]));
		}
		std::cout << n << " " << norm2 << " " << norm_max << std::endl;
	}

	//CHECK(norm2 == Approx(0.29078).margin(1e-6));
}

namespace{

class ACNConvDiffFemWorker{
	double nonstat_solution(Point p, double t) const {
		constexpr double pi = 3.141592653589793238462643;
		Vector vel = velocity(p);
		/*double r2 = (p.x() - vel.x()*t)*(p.x() - vel.x()*t);
		return 1.0/std::sqrt(4*pi*_eps*(t + _t0)) * exp(-r2/(4*_eps*(t + _t0)));*/

		////////////////////////////////////////////////////////////////////////////////////////////
		double x1 = 0.5 + (_x01 - 0.5) * cos(t) - (_y01 - 0.5) * sin(t);
		double y1 = 0.5 + (_x01 - 0.5) * sin(t) + (_y01 - 0.5) * cos(t);
		double r2 = (p.x() - x1) * (p.x() - x1) + (p.y() - y1) * (p.y() - y1);
		return 1.0 / (4 * pi * _eps * (t + _t0)) * exp(-r2 / (4 * _eps * (t + _t0)));
		////////////////////////////////////////////////////////////////////////////////////////////
	}
public:
	Vector velocity(Point p) const{
		/*return {1, 0, 0};*/

		////////////////////////////////////////////////////////////////////////////////////////////
		return { - p.y() + 0.5, p.x() - 0.5, 0};
		////////////////////////////////////////////////////////////////////////////////////////////
	};
	double exact_solution(Point p) const{
		return nonstat_solution(p, _time);
	}

	ACNConvDiffFemWorker(const IGrid& g, double eps, double tau)
			: _grid(g), _eps(eps), _tau(tau),
			  _fem(build_fem(g)), _u(_fem.n_bases(), 0){};
	
	void initialize(){
		// velocity
		for (size_t i=0; i<_fem.n_bases(); ++i){
			Vector v = velocity(_fem.reference_point(i));
			_vx.push_back(v.x());
			_vy.push_back(v.y());
			_vz.push_back(v.z());
		}

		// initial solution
		for (size_t i=0; i<_fem.n_bases(); ++i){
			_u[i] = exact_solution(_fem.reference_point(i));
		}

		// lhs, rhs
		assemble_solver();
	}

	void step(){
		_time += _tau;
		impl_step();
	}

	void save_vtk(const std::string& filename){
		// save grid
		_grid.save_vtk(filename);
		
		// save numerical solution
		VtkUtils::add_point_data(_u, "numerical", filename, _grid.n_points());

		// save exact solution
		std::vector<double> exact(_grid.n_points());
		for (size_t i=0; i<_grid.n_points(); ++i){
			exact[i] = exact_solution(_grid.point(i));
		}
		VtkUtils::add_point_data(exact, "exact", filename, _grid.n_points());
	}

	double current_time() const {
		return _time;
	}

	double compute_norm2() const{
		std::vector<double> diff2;
		double vol = 1.0;
		for (size_t i=0; i<_fem.n_bases(); ++i){
			double diff = _u[i] - exact_solution(_fem.reference_point(i));
			diff2.push_back(diff*diff);
		}
		diff2 = _M.mult_vec(diff2);
		double sum = std::accumulate(diff2.begin(), diff2.end(), 0.0);
		return std::sqrt(sum / vol);
	}

protected:
	const IGrid& _grid;
	const double _eps;
	const double _tau;
	FemAssembler _fem;
	std::vector<double> _u, _vx, _vy, _vz;
	CsrMatrix _M, _Rhs;
	double _time = 0.0;
	//double _t0 = 3.0;

	//////////////////////////////////////////////////////////////////////////////////////
	double _t0 = 0.3;
	//////////////////////////////////////////////////////////////////////////////////////

	AmgcMatrixSolver _solver;

	/////////////////////////////////////////////////////////////////////////////////////
	double _x01 = 0.8; double _y01 = 0.5;
	/////////////////////////////////////////////////////////////////////////////////////

	static FemAssembler build_fem(const IGrid& grid){
		size_t n_bases = grid.n_points();
		std::vector<FemElement> elements;
		std::vector<std::vector<size_t>> tab_elem_basis;

		// elements
		for (size_t icell=0; icell < grid.n_cells(); ++icell){
			std::vector<size_t> ipoints = grid.tab_cell_point(icell);
			tab_elem_basis.push_back(ipoints);
			if (ipoints.size() == 2){
				Point p0 = grid.point(ipoints[0]);
				Point p1 = grid.point(ipoints[1]);
				
				auto geom = std::make_shared<SegmentLinearGeometry>(p0, p1);
				auto basis = std::make_shared<SegmentLinearBasis>();
				const Quadrature* quadrature = quadrature_segment_gauss2();
				auto integrals = std::make_shared<NumericElementIntegrals>(quadrature, geom, basis);
				elements.push_back(FemElement{geom, basis, integrals});
			} else if (ipoints.size() == 3){
				Point p0 = grid.point(ipoints[0]);
				Point p1 = grid.point(ipoints[1]);
				Point p2 = grid.point(ipoints[2]);
				
				auto geom = std::make_shared<TriangleLinearGeometry>(p0, p1, p2);
				auto basis = std::make_shared<TriangleLinearBasis>();
				const Quadrature* quadrature = quadrature_triangle_gauss2();
				auto integrals = std::make_shared<NumericElementIntegrals>(quadrature, geom, basis);
				elements.push_back(FemElement{geom, basis, integrals});
			} else if (ipoints.size() == 4){
				Point p0 = grid.point(ipoints[0]);
				Point p1 = grid.point(ipoints[1]);
				Point p2 = grid.point(ipoints[2]);
				Point p3 = grid.point(ipoints[3]);
				
				auto geom = std::make_shared<QuadrangleLinearGeometry>(p0, p1, p2, p3);
				auto basis = std::make_shared<QuadrangleLinearBasis>();
				const Quadrature* quadrature = quadrature_square_gauss2();
				auto integrals = std::make_shared<NumericElementIntegrals>(quadrature, geom, basis);
				elements.push_back(FemElement{geom, basis, integrals});
			} else {
				throw std::runtime_error("invalid fem grid");
			}
		}

		return FemAssembler(n_bases, elements, tab_elem_basis);
	}

	void impl_step(){
		std::vector<double> rhs = assemble_rhs();
		_solver.solve(rhs, _u);
	}

	virtual void assemble_solver() = 0;

	std::vector<double> assemble_rhs(){
		// matrix mult
		std::vector<double> ret = _Rhs.mult_vec(_u);
		// boundary_condition
		/*ret[0] = exact_solution(_grid.point(0));
		ret[_grid.n_points()-1] = exact_solution(_grid.point(_grid.n_points()-1));*/

		/////////////////////////////////////////////////////////////////////////////////////////////////////
		std::vector<size_t> boundaries = _grid.boundary_points();
		for (size_t i : boundaries)
		{
			ret[i] = exact_solution(_grid.point(i));
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////

		return ret;
	}
};

class SupgWorker: public ACNConvDiffFemWorker{
public:
	SupgWorker(const IGrid& grid, double eps, double tau, double s_supg)
		:ACNConvDiffFemWorker(grid, eps, tau), _supg_coef(s_supg){}

protected:
	const double _supg_coef;

	void assemble_solver() override{
		// Galerkin
		CsrMatrix M(_fem.stencil());
		CsrMatrix K(_fem.stencil());
		CsrMatrix D(_fem.stencil());
		// SUPG
		CsrMatrix Ms(_fem.stencil());
		CsrMatrix Ks(_fem.stencil());
		CsrMatrix Ds(_fem.stencil());
		for (size_t ielem=0; ielem < _fem.n_elements(); ++ielem){
			const FemElement& elem = _fem.element(ielem);
			auto vx = _fem.local_vector(ielem, _vx);
			auto vy = _fem.local_vector(ielem, _vy);
			auto vz = _fem.local_vector(ielem, _vz);
			// ====================== Galerkin
			// mass
			{
				std::vector<double> local = elem.integrals->mass_matrix();
				_fem.add_to_global_matrix(ielem, local, M.vals());
			}
			// transport
			{
				std::vector<double> local = elem.integrals->transport_matrix(vx, vy, vz);
				_fem.add_to_global_matrix(ielem, local, K.vals());
			}
			// diffusion
			{
				std::vector<double> local = elem.integrals->stiff_matrix();
				_fem.add_to_global_matrix(ielem, local, D.vals());
			}
			// ======================== Stabilization
			// mass
			{
				std::vector<double> local = elem.integrals->mass_matrix_stab_supg(vx, vy, vz);
				_fem.add_to_global_matrix(ielem, local, Ms.vals());
			}
			// transport
			{
				std::vector<double> local = elem.integrals->transport_matrix_stab_supg(vx, vy, vz);
				_fem.add_to_global_matrix(ielem, local, Ks.vals());
			}
			// diffusion
			{
				std::vector<double> local = elem.integrals->stiff_matrix_stab_supg(vx, vy, vz);
				_fem.add_to_global_matrix(ielem, local, Ds.vals());
			}
		}
		// slae matrix
		double theta = 0.5;
		CsrMatrix Lhs(_fem.stencil());
		CsrMatrix Rhs(_fem.stencil());
		for (size_t a=0; a<Lhs.n_nonzeros(); ++a){
			Lhs.vals()[a] =
				M.vals()[a]
				+ _tau*theta*K.vals()[a]
				+ _tau*theta*_eps * D.vals()[a]

				+ _supg_coef*Ms.vals()[a]
				+ _supg_coef*_tau*theta*Ks.vals()[a]
				+ _supg_coef*_tau*theta*_eps*Ds.vals()[a]
				;
			Rhs.vals()[a] =
				M.vals()[a]
				- _tau*(1 - theta)*K.vals()[a]
				- _tau*(1 - theta)*_eps*D.vals()[a]

				+ _supg_coef*Ms.vals()[a]
				- _supg_coef*_tau*(1 - theta)*Ks.vals()[a]
				- _supg_coef*_tau*(1 - theta)*_eps*Ds.vals()[a]
				;
		}

		// Boundary conditions
		//Lhs.set_unit_row(0);
		//Lhs.set_unit_row(_grid.n_points()-1);

		////////////////////////////////////////////////////////////////////////////////////////////////
		std::vector<size_t> boundaries = _grid.boundary_points();
		for (size_t i : boundaries)
		{
			Lhs.set_unit_row(i);
		}
		////////////////////////////////////////////////////////////////////////////////////////////////

		// Solver
		_solver.set_matrix(Lhs);

		_M = std::move(M);
		_Rhs = std::move(Rhs);
	}
};

}

TEST_CASE("1D convection-diffusion with SUPG", "[convdiff-fem-supg]"){
	std::cout << std::endl << "--- cfd24_test [convdiff-fem-supg] --- " << std::endl;
	double tend = 2.0;
	double h = 0.1;
	double Lx = 4;
	double Cu = 0.5;
	double eps = 1e-3;
	double s_supg = 1e-2;

	// solver
	Grid1D grid(0, Lx, Lx / h);
	double tau = Cu * h;
	SupgWorker worker(grid, eps, tau, s_supg);
	worker.initialize();

	// saver
	VtkUtils::TimeSeriesWriter writer("convdiff-supg");
	std::string out_filename = writer.add(worker.current_time());
	worker.save_vtk(out_filename);

	double n2;
	while (worker.current_time() < tend - 1e-6) {
		// solve problem
		worker.step();
		// export solution to vtk
		out_filename = writer.add(worker.current_time());
		worker.save_vtk(out_filename);

		n2 = worker.compute_norm2();

		std::cout << worker.current_time() << " " << n2 << std::endl;
	};
	CHECK(n2 == Approx(0.615742).margin(1e-6));
}

namespace {

class CgWorker: public ACNConvDiffFemWorker{
public:
	CgWorker(const IGrid& grid, double eps, double tau)
		:ACNConvDiffFemWorker(grid, eps, tau){}

protected:
	void assemble_solver() override{
		// Galerkin
		CsrMatrix M(_fem.stencil());
		CsrMatrix K(_fem.stencil());
		CsrMatrix D(_fem.stencil());
		// CG
		CsrMatrix Ks(_fem.stencil());
		for (size_t ielem=0; ielem < _fem.n_elements(); ++ielem){
			const FemElement& elem = _fem.element(ielem);
			auto vx = _fem.local_vector(ielem, _vx);
			auto vy = _fem.local_vector(ielem, _vy);
			auto vz = _fem.local_vector(ielem, _vz);
			// ====================== Galerkin
			// mass
			{
				std::vector<double> local = elem.integrals->mass_matrix();
				_fem.add_to_global_matrix(ielem, local, M.vals());
			}
			// transport
			{
				std::vector<double> local = elem.integrals->transport_matrix(vx, vy, vz);
				_fem.add_to_global_matrix(ielem, local, K.vals());
			}
			// diffusion
			{
				std::vector<double> local = elem.integrals->stiff_matrix();
				_fem.add_to_global_matrix(ielem, local, D.vals());
			}
			// ======================== Stabilization
			// transport
			{
				std::vector<double> local = elem.integrals->transport_matrix_stab_supg(vx, vy, vz);
				_fem.add_to_global_matrix(ielem, local, Ks.vals());
			}
		}
		// slae matrix
		double theta = 0.5;
		CsrMatrix Lhs(_fem.stencil());
		CsrMatrix Rhs(_fem.stencil());
		for (size_t a=0; a<Lhs.n_nonzeros(); ++a){
			Lhs.vals()[a] =
				M.vals()[a]
				+ _tau*theta*_eps * D.vals()[a]
				;
			Rhs.vals()[a] =
				M.vals()[a]
				- _tau*K.vals()[a]
				- _tau*(1 - theta)*_eps*D.vals()[a]
				- _tau*_tau/2.0*Ks.vals()[a]
				;
		}

		// Boundary conditions
		//Lhs.set_unit_row(0);
		//Lhs.set_unit_row(_grid.n_points()-1);

		std::vector<size_t> boundaries = _grid.boundary_points();
		for (size_t i : boundaries)
		{
			Lhs.set_unit_row(i);
		}
		
		// Solver
		_solver.set_matrix(Lhs);

		_M = std::move(M);
		_Rhs = std::move(Rhs);
	}
};
}

TEST_CASE("1D convection-diffusion with CG", "[convdiff-fem-cg]"){
	std::cout << std::endl << "--- cfd24_test [convdiff-fem-cg] --- " << std::endl;
	double tend = 2.0;
	double h = 0.1;
	double Lx = 4;
	double Cu = 0.1;
	double eps = 1e-3;

	// solver
	Grid1D grid(0, Lx, Lx / h);
	double tau = Cu * h;
	/////////////////////////////////////////////////////////////////////////////////////
	//std::cout << "tau = " << tau << std::endl;
	/////////////////////////////////////////////////////////////////////////////////////
	CgWorker worker(grid, eps, tau);
	worker.initialize();

	// saver
	VtkUtils::TimeSeriesWriter writer("convdiff-cg");
	std::string out_filename = writer.add(worker.current_time());
	worker.save_vtk(out_filename);

	double n2;
	while (worker.current_time() < tend - 1e-6) {
		// solve problem
		worker.step();
		// export solution to vtk
		out_filename = writer.add(worker.current_time());
		worker.save_vtk(out_filename);

		n2 = worker.compute_norm2();

		std::cout << worker.current_time() << " " << n2 << std::endl;
	};
	CHECK(n2 == Approx(0.937086).margin(1e-6));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("2D convection-diffusion with CG", "[convdiff-fem-cg-2D]") {
	std::cout << std::endl << "--- cfd24_test [convdiff-fem-cg-2D] --- " << std::endl;
	double tend = 6.283185307179;
	double h = 1.0/sqrt(5000);
	double Lx = 1.0;
	double Cu = 0.5;
	double eps = 1e-3;

	// solver
	//RegularGrid2D grid(0.0, 1.0, 0.0, 1.0, 1.0/h, 1.0/h);
	//double tau = Cu * h;
	double tau = 0.005;
	/////////////////////////////////////////////////////////////////////////////////////
	//std::cout << "tau = " << tau << std::endl;
	
	std::string fn = test_directory_file("tetragrid.vtk");
	auto grid = UnstructuredGrid2D::vtk_read(fn);
	 
	/////////////////////////////////////////////////////////////////////////////////////
	CgWorker worker(grid, eps, tau);
	worker.initialize();

	// saver
	VtkUtils::TimeSeriesWriter writer("convdiff-cg-tetra");
	std::string out_filename = writer.add(worker.current_time());
	worker.save_vtk(out_filename);

	double n2;
	while (worker.current_time() < tend - 1e-6) {
		// solve problem
		worker.step();
		// export solution to vtk
		out_filename = writer.add(worker.current_time());
		worker.save_vtk(out_filename);

		n2 = worker.compute_norm2();

		std::cout << worker.current_time() << " " << n2 << std::endl;
	};
	CHECK(n2 == Approx(0.937086).margin(1e-6));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////