#include "cfd24_test.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/mat/lodmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/debug/printer.hpp"
#include "utils/vecmat.hpp"
#include "cfd24/debug/tictoc.hpp"
#include "cfd24/debug/saver.hpp"
#include "cfd24/fvm/fvm_assembler.hpp"
#include "utils/filesystem.hpp"
#include "cfd24/grid/unstructured_grid2d.hpp"
#include "cfd24/fvm/fvm_dpdn_boundary.hpp"
#include <list>

using namespace cfd;

struct Cavern2DFvmSimpleWorker {
	Cavern2DFvmSimpleWorker(const IGrid& grid, double Re, double E);
	void initialize_saver(std::string stem);
	double step();
	void save_current_fields(size_t iter);

	size_t vec_size() const {
		return _collocations.size();
	}
private:
	const IGrid& _grid;
	const double _Re;
	const double _tau;
	const double _alpha_p;
	const FvmExtendedCollocations _collocations;
	const FvmFacesDn _dfdn_computer;
	const FvmCellGradient _grad_computer;
	FvmDpDnComputer2D _dpdn_computer;

	struct BoundaryInfo {
		std::vector<size_t> bnd_colloc;
		std::vector<size_t> bnd_colloc_u0;
		std::vector<size_t> bnd_colloc_u1;
	};
	BoundaryInfo _boundary_info;

	std::vector<double> _p;
	std::vector<double> _u;
	std::vector<double> _v;
	std::vector<double> _un_face;
	std::vector<double> _dpdn_face;
	std::vector<Vector> _grad_p;

	AmgcMatrixSolver _p_stroke_solver;
	AmgcMatrixSolver _uv_solver;

	CsrMatrix _mat_uv;
	std::vector<double> _rhs_u;
	std::vector<double> _rhs_v;
	double _d;

	std::shared_ptr<VtkUtils::TimeSeriesWriter> _writer;

	double to_next_iteration();
	void gather_boundary_collocations();
	void assemble_p_stroke_solver();
	void assemble_uv_slae();

	std::vector<double> compute_u_star();
	std::vector<double> compute_v_star();
	std::vector<double> compute_un_star_face_rhie_chow(const std::vector<double>& u_star, const std::vector<double>& v_star);
	std::vector<double> compute_un_stroke_face(const std::vector<double>& dpstroke_dn);
	std::vector<double> compute_p_stroke(const std::vector<double>& un_star_face);
	std::vector<double> compute_u_stroke(const std::vector<Vector>& grad_p_stroke);
	std::vector<double> compute_v_stroke(const std::vector<Vector>& grad_p_stroke);

	static double compute_tau(const IGrid& grid, double Re, double E);
};

double Cavern2DFvmSimpleWorker::compute_tau(const IGrid& grid, double Re, double E) {
	double h2 = grid.cell_volume(0);
	for (size_t i = 1; i < grid.n_cells(); ++i) {
		h2 = std::max(h2, grid.cell_volume(i));
	}
	return E * Re * h2 / 4.0;
}

Cavern2DFvmSimpleWorker::Cavern2DFvmSimpleWorker(const IGrid& grid, double Re, double E) :
	_grid(grid),
	_Re(Re),
	_tau(compute_tau(grid, Re, E)),
	_alpha_p(1.0 / (1.0 + E)),
	_collocations(grid),
	_dfdn_computer(grid, _collocations),
	_grad_computer(grid, _collocations),
	_dpdn_computer(Re, _alpha_p, _tau, 1.0 / (1 + E), grid, _collocations)
{
	_d = 1.0 / (1 + E);
	gather_boundary_collocations();
	assemble_p_stroke_solver();

	_u = std::vector<double>(vec_size(), 0);
	_v = std::vector<double>(vec_size(), 0);
	_p = std::vector<double>(vec_size(), 0);
	_dpdn_face = std::vector<double>(_grid.n_faces(), 0);
	_un_face = std::vector<double>(_grid.n_faces(), 0);
	_grad_p = std::vector<Vector>(_grid.n_cells(), { 0, 0, 0 });
	to_next_iteration();
}

void Cavern2DFvmSimpleWorker::gather_boundary_collocations() {
	BoundaryInfo& bi = _boundary_info;

	std::list<std::pair<size_t, size_t>> colloc_faces;
	for (size_t icolloc : _collocations.face_collocations) {
		size_t iface = _collocations.face_index(icolloc);
		bi.bnd_colloc.push_back(icolloc);
		if (std::abs(_grid.face_center(iface).y() - 1) < 1e-6) {
			bi.bnd_colloc_u1.push_back(icolloc);
		}
		else {
			bi.bnd_colloc_u0.push_back(icolloc);
		}
	}
}

void Cavern2DFvmSimpleWorker::initialize_saver(std::string stem) {
	_writer.reset(new VtkUtils::TimeSeriesWriter(stem));
};

double Cavern2DFvmSimpleWorker::to_next_iteration() {
	assemble_uv_slae();

	// residual vectors
	std::vector<double> res_u = compute_residual_vec(_mat_uv, _rhs_u, _u);
	std::vector<double> res_v = compute_residual_vec(_mat_uv, _rhs_v, _v);
	for (size_t icell = 0; icell < _grid.n_cells(); ++icell) {
		double coef = 1.0 / _tau / _grid.cell_volume(icell);
		res_u[icell] *= coef;
		res_v[icell] *= coef;
	}

	// norm
	double res = 0;
	for (size_t icell = 0; icell < _grid.n_cells(); ++icell) {
		res = std::max(res, std::max(res_u[icell], res_v[icell]));
	}
	return res;
};

double Cavern2DFvmSimpleWorker::step() {
	// Predictor step: U-star
	std::vector<double> u_star = compute_u_star();
	std::vector<double> v_star = compute_v_star();
	std::vector<double> un_star_face = compute_un_star_face_rhie_chow(u_star, v_star);
	// Pressure correction
	std::vector<double> p_stroke = compute_p_stroke(un_star_face);
	std::vector<Vector> grad_p_stroke = _grad_computer.compute(p_stroke);
	std::vector<double> dpstroke_dn_face = _dfdn_computer.compute(p_stroke);
	// Velocity correction
	std::vector<double> u_stroke = compute_u_stroke(grad_p_stroke);
	std::vector<double> v_stroke = compute_v_stroke(grad_p_stroke);
	std::vector<double> un_stroke_face = compute_un_stroke_face(dpstroke_dn_face);

	// Set final values
	_u = vector_sum(u_star, 1.0, u_stroke);
	_v = vector_sum(v_star, 1.0, v_stroke);
	_un_face = vector_sum(un_star_face, 1.0, un_stroke_face);
	_p = vector_sum(_p, _alpha_p, p_stroke);
	_grad_p = vector_sum(_grad_p, _alpha_p, grad_p_stroke);
	_dpdn_face = vector_sum(_dpdn_face, _alpha_p, dpstroke_dn_face);

	return to_next_iteration();
}

void Cavern2DFvmSimpleWorker::save_current_fields(size_t iter) {
	if (_writer) {
		std::string filepath = _writer->add(iter);
		_grid.save_vtk(filepath);
		VtkUtils::add_cell_data(_p, "pressure", filepath, _grid.n_cells());
		VtkUtils::add_cell_vector(_u, _v, "velocity", filepath, _grid.n_cells());
	}
}

void Cavern2DFvmSimpleWorker::assemble_p_stroke_solver() {
	LodMatrix mat(vec_size());
	for (size_t iface = 0; iface < _grid.n_faces(); ++iface) {
		size_t negative_colloc = _collocations.tab_face_colloc(iface)[0];
		size_t positive_colloc = _collocations.tab_face_colloc(iface)[1];
		double area = _grid.face_area(iface);

		for (const std::pair<const size_t, double>& iter : _dfdn_computer.linear_combination(iface)) {
			size_t column = iter.first;
			double coef = -_d * area * iter.second;
			mat.add_value(negative_colloc, column, coef);
			mat.add_value(positive_colloc, column, -coef);
		}
	}
	mat.set_unit_row(0);
	_p_stroke_solver.set_matrix(mat.to_csr());
}

void Cavern2DFvmSimpleWorker::assemble_uv_slae() {
	// =============== Left side
	LodMatrix mat(vec_size());
	// u
	for (size_t icell = 0; icell < _grid.n_cells(); ++icell) {
		mat.add_value(icell, icell, _grid.cell_volume(icell));
	}
	for (size_t iface = 0; iface < _grid.n_faces(); ++iface) {
		size_t negative_colloc = _collocations.tab_face_colloc(iface)[0];
		size_t positive_colloc = _collocations.tab_face_colloc(iface)[1];
		double area = _grid.face_area(iface);

		// - tau/Re * Laplace(u)
		for (const std::pair<const size_t, double>& iter : _dfdn_computer.linear_combination(iface)) {
			size_t column = iter.first;
			double coef = -_tau / _Re * area * iter.second;
			mat.add_value(negative_colloc, column, coef);
			mat.add_value(positive_colloc, column, -coef);
		}

		// + convection
		{
			double coef = _tau * area * _un_face[iface] / 2.0;
			mat.add_value(negative_colloc, negative_colloc, coef);
			mat.add_value(negative_colloc, positive_colloc, coef);
			mat.add_value(positive_colloc, positive_colloc, -coef);
			mat.add_value(positive_colloc, negative_colloc, -coef);
		}
	}
	// bnd
	for (size_t icolloc : _boundary_info.bnd_colloc) {
		mat.set_unit_row(icolloc);
	}
	_mat_uv = mat.to_csr();
	_uv_solver.set_matrix(_mat_uv);

	// ============== right side
	_rhs_u.resize(vec_size());
	_rhs_v.resize(vec_size());
	for (size_t icell = 0; icell < _grid.n_cells(); ++icell) {
		_rhs_u[icell] = (_u[icell] - _tau * _grad_p[icell].x()) * _grid.cell_volume(icell);
		_rhs_v[icell] = (_v[icell] - _tau * _grad_p[icell].y()) * _grid.cell_volume(icell);
	}
	// bnd
	for (size_t icolloc : _boundary_info.bnd_colloc_u0) {
		_rhs_u[icolloc] = 0;
		_rhs_v[icolloc] = 0;
	}
	for (size_t icolloc : _boundary_info.bnd_colloc_u1) {
		_rhs_u[icolloc] = 1;
		_rhs_v[icolloc] = 0;
	}
}

std::vector<double> Cavern2DFvmSimpleWorker::compute_u_star() {
	std::vector<double> u_star(_u);
	_uv_solver.solve(_rhs_u, u_star);
	return u_star;
}

std::vector<double> Cavern2DFvmSimpleWorker::compute_v_star() {
	std::vector<double> v_star(_v);
	_uv_solver.solve(_rhs_v, v_star);
	return v_star;
}

std::vector<double> Cavern2DFvmSimpleWorker::compute_un_star_face_rhie_chow(
	const std::vector<double>& u_star,
	const std::vector<double>& v_star) {

	std::vector<double> ret(_grid.n_faces());

	for (size_t iface = 0; iface < _grid.n_faces(); ++iface) {
		size_t ci = _grid.tab_face_cell(iface)[0];
		size_t cj = _grid.tab_face_cell(iface)[1];
		if (ci == INVALID_INDEX || cj == INVALID_INDEX) {
			ret[iface] = 0;
		}
		else {
			// Rhie-Chow interpolation
			Vector normal = _grid.face_normal(iface);
			Vector uvec_i = Vector(u_star[ci], v_star[ci]);
			Vector uvec_j = Vector(u_star[cj], v_star[cj]);
			//Vector uvec_old_i = Vector(_u[ci], _v[ci]);
			//Vector uvec_old_j = Vector(_u[cj], _v[cj]);

			double ustar_i = dot_product(uvec_i, normal);
			double ustar_j = dot_product(uvec_j, normal);
			double dpdn_i = dot_product(_grad_p[ci], normal);
			double dpdn_j = dot_product(_grad_p[cj], normal);
			double dpdn_ij = _dpdn_face[iface];
			//double u_i = dot_product(uvec_old_i, normal);
			//double u_j = dot_product(uvec_old_j, normal);
			//double u_ij = _un_face[iface];

			ret[iface] = 0.5 * (ustar_i + ustar_j)
				+ 0.5 * _tau * (dpdn_i + dpdn_j - 2 * dpdn_ij);
			//+ 0.5*(2*u_ij - u_i - u_j);
		}
	}

	return ret;
}

std::vector<double> Cavern2DFvmSimpleWorker::compute_un_stroke_face(const std::vector<double>& dpstroke_dn) {
	std::vector<double> un(_grid.n_faces());
	for (size_t iface = 0; iface < _grid.n_faces(); ++iface) {
		un[iface] = -_tau * _d * dpstroke_dn[iface];
	}
	return un;
}

std::vector<double> Cavern2DFvmSimpleWorker::compute_p_stroke(const std::vector<double>& un_star_face) {
	std::vector<double> rhs(vec_size(), 0.0);
	// internal
	for (size_t iface = 0; iface < _grid.n_faces(); ++iface) {
		double coef = -_grid.face_area(iface) / _tau * un_star_face[iface];
		size_t neg = _grid.tab_face_cell(iface)[0];
		size_t pos = _grid.tab_face_cell(iface)[1];
		if (neg != INVALID_INDEX) rhs[neg] += coef;
		if (pos != INVALID_INDEX) rhs[pos] -= coef;
	}
	rhs[0] = 0;
	// solve
	std::vector<double> p_stroke;
	_p_stroke_solver.solve(rhs, p_stroke);
	return p_stroke;
}

std::vector<double> Cavern2DFvmSimpleWorker::compute_u_stroke(const std::vector<Vector>& grad_p_stroke) {
	std::vector<double> u_stroke(vec_size(), 0.0);
	for (size_t i = 0; i < _grid.n_cells(); ++i) {
		u_stroke[i] = -_tau * _d * grad_p_stroke[i].x();
	}
	return u_stroke;
}

std::vector<double> Cavern2DFvmSimpleWorker::compute_v_stroke(const std::vector<Vector>& grad_p_stroke) {
	std::vector<double> v_stroke(vec_size(), 0.0);
	for (size_t i = 0; i < _grid.n_cells(); ++i) {
		v_stroke[i] = -_tau * _d * grad_p_stroke[i].y();
	}
	return v_stroke;
}

TEST_CASE("Cavern 2D, FVM-SIMPLE algorithm", "[cavern2-fvm-simple]") {
	std::cout << std::endl << "--- cfd24_test [cavern2-fvm-simple] --- " << std::endl;

	// problem parameters
	double Re = 100;
	size_t max_it = 10'000;
	double eps = 1e-0;
	double E = 2;

	// worker initialization
	RegularGrid2D grid(0, 1, 0, 1, 30, 30);
	Cavern2DFvmSimpleWorker worker(grid, Re, E);
	worker.initialize_saver("cavern2-fvm");

	// iterations loop
	size_t it = 0;
	for (it = 1; it < max_it; ++it) {
		double nrm = worker.step();

		// print norm and pressure value at the top-right corner
		std::cout << it << " " << nrm << std::endl;

		// export solution to vtk
		worker.save_current_fields(it);

		// break if residual is low enough
		if (nrm < eps) {
			break;
		}
	}
	CHECK(it == 6);
}