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

struct Cavern2DFemSimpleWorker{
	Cavern2DFemSimpleWorker(const IGrid& grid, double Re, double E);
	void initialize_saver(std::string stem);
	double step();
	void save_current_fields(size_t iter);
private:
	const IGrid& _grid;
	const double _Re;
	const double _tau;
	const double _alpha_p;

	std::vector<double> _p;
	std::vector<double> _u;
	std::vector<double> _v;

	AmgcMatrixSolver _p_stroke_solver;
	AmgcMatrixSolver _uv_solver;

	CsrMatrix _mat_uv;
	std::vector<double> _rhs_u;
	std::vector<double> _rhs_v;
	double _d;

	std::shared_ptr<VtkUtils::TimeSeriesWriter> _writer;

	size_t u_size() const;
	size_t p_size() const;

	double to_next_iteration();
	void assemble_p_stroke_solver();
	void assemble_uv_slae();

	std::vector<double> compute_u_star();
	std::vector<double> compute_v_star();
	std::vector<double> compute_p_stroke(const std::vector<double>& u_star, const std::vector<double>& v_star);
	std::vector<double> compute_u_stroke(const std::vector<double>& p_stroke);
	std::vector<double> compute_v_stroke(const std::vector<double>& p_stroke);

	static double compute_tau(const IGrid& grid, double Re, double E);
};

size_t Cavern2DFemSimpleWorker::u_size() const{
	_THROW_NOT_IMP_;
}

size_t Cavern2DFemSimpleWorker::p_size() const{
	_THROW_NOT_IMP_;
}

double Cavern2DFemSimpleWorker::compute_tau(const IGrid& grid, double Re, double E){
	double h2 = grid.cell_volume(0);
	for (size_t i=1; i<grid.n_cells(); ++i){
		h2 = std::max(h2, grid.cell_volume(i));
	}
	return E*Re*h2/4.0;
}

Cavern2DFemSimpleWorker::Cavern2DFemSimpleWorker(const IGrid& grid, double Re, double E):
	_grid(grid),
	_Re(Re),
	_tau(compute_tau(grid, Re, E)),
	_alpha_p(1.0/(1.0 + E))
{
	_d = 1.0/(1 + E);
	assemble_p_stroke_solver();

	_u = std::vector<double>(u_size(), 0);
	_v = std::vector<double>(u_size(), 0);
	_p = std::vector<double>(p_size(), 0);
	to_next_iteration();
}

void Cavern2DFemSimpleWorker::initialize_saver(std::string stem){
	_writer.reset(new VtkUtils::TimeSeriesWriter(stem));
};

double Cavern2DFemSimpleWorker::to_next_iteration(){
	assemble_uv_slae();

	// residual vectors
	std::vector<double> res_u = compute_residual_vec(_mat_uv, _rhs_u, _u);
	std::vector<double> res_v = compute_residual_vec(_mat_uv, _rhs_v, _v);
	for (size_t icell=0; icell < _grid.n_cells(); ++icell){
		double coef = 1.0 / _tau / _grid.cell_volume(icell);
		res_u[icell] *= coef;
		res_v[icell] *= coef;
	}

	// norm
	double res = 0;
	for (size_t icell=0; icell < _grid.n_cells(); ++icell){
		res = std::max(res, std::max(res_u[icell], res_v[icell]));
	}
	return res;
};

double Cavern2DFemSimpleWorker::step(){
	// Predictor step: U-star
	std::vector<double> u_star = compute_u_star();
	std::vector<double> v_star = compute_v_star();
	// Pressure correction
	std::vector<double> p_stroke = compute_p_stroke(u_star, v_star);
	// Velocity correction
	std::vector<double> u_stroke = compute_u_stroke(p_stroke);
	std::vector<double> v_stroke = compute_v_stroke(p_stroke);

	// Set final values
	_u = vector_sum(u_star, 1.0, u_stroke);
	_v= vector_sum(v_star, 1.0, v_stroke);
	_p = vector_sum(_p, _alpha_p, p_stroke);

	return to_next_iteration();
}

void Cavern2DFemSimpleWorker::save_current_fields(size_t iter){
	if (_writer){
		std::string filepath = _writer->add(iter);
		_grid.save_vtk(filepath);
		VtkUtils::add_cell_data(_p, "pressure", filepath, _grid.n_cells());
		VtkUtils::add_cell_vector(_u, _v, "velocity", filepath, _grid.n_cells());
	}
}

void Cavern2DFemSimpleWorker::assemble_p_stroke_solver(){
	LodMatrix mat(p_size());
	_THROW_NOT_IMP_;
	mat.set_unit_row(0);
	_p_stroke_solver.set_matrix(mat.to_csr());
}

void Cavern2DFemSimpleWorker::assemble_uv_slae(){
	// =============== Left side
	LodMatrix mat(u_size());
	_THROW_NOT_IMP_;
	_mat_uv = mat.to_csr();
	_uv_solver.set_matrix(_mat_uv);

	// ============== right side
	_rhs_u.resize(u_size());
	_rhs_v.resize(u_size());
}

std::vector<double> Cavern2DFemSimpleWorker::compute_u_star(){
	std::vector<double> u_star(_u);
	_uv_solver.solve(_rhs_u, u_star);
	return u_star;
}

std::vector<double> Cavern2DFemSimpleWorker::compute_v_star(){
	std::vector<double> v_star(_v);
	_uv_solver.solve(_rhs_v, v_star);
	return v_star;
}

std::vector<double> Cavern2DFemSimpleWorker::compute_p_stroke(const std::vector<double>& u_star, const std::vector<double>& v_star){
	std::vector<double> rhs(p_size(), 0.0);
	_THROW_NOT_IMP_;
	// solve
	std::vector<double> p_stroke;
	_p_stroke_solver.solve(rhs, p_stroke);
	return p_stroke;
}

std::vector<double> Cavern2DFemSimpleWorker::compute_u_stroke(const std::vector<double>& p_stroke){
	std::vector<double> u_stroke(u_size(), 0.0);
	_THROW_NOT_IMP_;
	return u_stroke;
}

std::vector<double> Cavern2DFemSimpleWorker::compute_v_stroke(const std::vector<double>& p_stroke){
	std::vector<double> v_stroke(u_size(), 0.0);
	_THROW_NOT_IMP_;
	return v_stroke;
}

TEST_CASE("Cavern 2D, FEM-SIMPLE algorithm", "[cavern2-fem-simple]"){
	std::cout << std::endl << "--- cfd24_test [cavern2-fem-simple] --- " << std::endl;

	// problem parameters
	double Re = 100;
	size_t max_it = 10'000;
	double eps = 1e-0;
	double E = 2;

	// worker initialization
	RegularGrid2D grid(0, 1, 0, 1, 30, 30);
	Cavern2DFemSimpleWorker worker(grid, Re, E);
	worker.initialize_saver("cavern2-fem");

	// iterations loop
	size_t it = 0;
	for (it=1; it < max_it; ++it){
		double nrm = worker.step();

		// print norm and pressure value at the top-right corner
		std::cout << it << " " << nrm << std::endl;

		// export solution to vtk
		worker.save_current_fields(it);

		// break if residual is low enough
		if (nrm < eps){
			break;
		}
	}
}
