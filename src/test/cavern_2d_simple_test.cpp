#include "cfd24_test.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/mat/lodmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/debug/printer.hpp"
#include "utils/vecmat.hpp"

#include <fstream>

using namespace cfd;

std::ofstream norms;

struct Cavern2DSimpleWorker{
	Cavern2DSimpleWorker(double Re, size_t n_cells, double tau, double alpha_p, double u_bottom);
	void initialize_saver(bool save_exact_fields, std::string stem);
	double set_uvp(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& p);
	double step();
	void save_current_fields(size_t iter);

	size_t u_size() const{
		return _yf_grid.n_points();
	}
	size_t v_size() const{
		return _xf_grid.n_points();
	}
	size_t p_size() const{
		return _cc_grid.n_points();
	}
	const std::vector<double>& pressure() const{
		return _p;
	}
private:
	const RegularGrid2D _grid;
	const RegularGrid2D _cc_grid;
	const RegularGrid2D _xf_grid;
	const RegularGrid2D _yf_grid;
	const double _hx;
	const double _hy;
	const double _Re;
	const double _tau;
	const double _alpha_p;

	std::vector<double> _p;
	std::vector<double> _u;
	std::vector<double> _v;

	double _du;
	double _dv;
	AmgcMatrixSolver _p_stroke_solver;

	CsrMatrix _mat_u;
	CsrMatrix _mat_v;
	std::vector<double> _rhs_u;
	std::vector<double> _rhs_v;

	std::shared_ptr<VtkUtils::TimeDependentWriter> _writer_u;
	std::shared_ptr<VtkUtils::TimeDependentWriter> _writer_v;
	std::shared_ptr<VtkUtils::TimeDependentWriter> _writer_p;
	std::shared_ptr<VtkUtils::TimeDependentWriter> _writer_all;

	void assemble_p_stroke_solver();
	void assemble_u_slae();
	void assemble_v_slae();

	std::vector<double> compute_u_star();
	std::vector<double> compute_v_star();
	std::vector<double> compute_p_stroke(const std::vector<double>& u_star, const std::vector<double>& v_star);
	std::vector<double> compute_u_stroke(const std::vector<double>& p_stroke);
	std::vector<double> compute_v_stroke(const std::vector<double>& p_stroke);
	double u_i_j(size_t i, size_t j) const;
	double u_ip_jp(size_t i, size_t j) const;
	double v_i_j(size_t i, size_t j) const;
	double v_ip_jp(size_t i, size_t j) const;
	double p_ip_jp(size_t i, size_t j) const;
	std::vector<Vector> build_main_grid_velocity() const;

	////////////////////////////////////////////////////////////////////
	double _u_bottom;
	virtual double get_bottom_velocity() const { return _u_bottom; };
	////////////////////////////////////////////////////////////////////
};

Cavern2DSimpleWorker::Cavern2DSimpleWorker(double Re, size_t n_cells, double tau, double alpha_p, double u_bottom):
	_grid(0, 1, 0, 1, n_cells, n_cells),
	_cc_grid(_grid.cell_centered_grid()),
	_xf_grid(_grid.xface_centered_grid()),
	_yf_grid(_grid.yface_centered_grid()),
	_hx(1.0/n_cells),
	_hy(1.0/n_cells),
	_Re(Re),
	_tau(tau),
	_alpha_p(alpha_p),
	_u_bottom(u_bottom)
{
	_du = 1.0 / (1 + 2.0*_tau/_Re * (1.0/_hx/_hx + 1.0/_hy/_hy));
	_dv = 1.0 / (1 + 2.0*_tau/_Re * (1.0/_hx/_hx + 1.0/_hy/_hy));
	assemble_p_stroke_solver();
}

void Cavern2DSimpleWorker::initialize_saver(bool save_exact_fields, std::string stem){
	_writer_all.reset(new VtkUtils::TimeDependentWriter(stem));
	if (save_exact_fields){
		_writer_u.reset(new VtkUtils::TimeDependentWriter(stem + "-u"));
		_writer_v.reset(new VtkUtils::TimeDependentWriter(stem + "-v"));
		_writer_p.reset(new VtkUtils::TimeDependentWriter(stem + "-p"));
	}
};

double Cavern2DSimpleWorker::set_uvp(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& p){
	_u = u;
	_v = v;
	_p = p;
	assemble_u_slae();
	assemble_v_slae();
	// residuals
	double nrm_u = compute_residual(_mat_u, _rhs_u, _u)/_tau;
	double nrm_v = compute_residual(_mat_v, _rhs_v, _v)/_tau;

	//for printing ru, rv
	//std::vector<double> ru = compute_residual2(_mat_u, _rhs_u, _u); 
	//std::vector<double> rv = compute_residual2(_mat_v, _rhs_v, _v);
	//
	//for (size_t i = 0; i < ru.size(); ++i)
	//{
	//	norms << ru[i] << "		" << rv[i] << std::endl;
	//}
	//norms << std::endl;
	
	return std::max(nrm_u, nrm_v);
};

double Cavern2DSimpleWorker::step(){
	// Predictor step: U-star
	std::vector<double> u_star = compute_u_star();
	std::vector<double> v_star = compute_v_star();
	// Pressure correction
	std::vector<double> p_stroke = compute_p_stroke(u_star, v_star);
	// Velocity correction
	std::vector<double> u_stroke = compute_u_stroke(p_stroke);
	std::vector<double> v_stroke = compute_v_stroke(p_stroke);
	// Set final values
	std::vector<double> u_new = vector_sum(u_star, 1.0, u_stroke);
	std::vector<double> v_new = vector_sum(v_star, 1.0, v_stroke);
	std::vector<double> p_new = vector_sum(_p, _alpha_p, p_stroke);

	return set_uvp(u_new, v_new, p_new);
}

void Cavern2DSimpleWorker::save_current_fields(size_t iter){
	// all data on the main grid
	if (_writer_all){
		std::string filepath = _writer_all->add(iter);
		_grid.save_vtk(filepath);
		VtkUtils::add_cell_data(_p, "pressure", filepath);
		VtkUtils::add_point_vector(build_main_grid_velocity(), "velocity", filepath);
	}
	// pressure
	if (_writer_p){
		std::string filepath = _writer_p->add(iter);
		_cc_grid.save_vtk(filepath);
		VtkUtils::add_point_data(_p, "pressure", filepath);
	}
	// u
	if (_writer_u){
		std::string filepath = _writer_u->add(iter);
		_yf_grid.save_vtk(filepath);
		VtkUtils::add_point_data(_u, "velocity-x", filepath);
	}
	// v
	if (_writer_v){
		std::string filepath = _writer_v->add(iter);
		_xf_grid.save_vtk(filepath);
		VtkUtils::add_point_data(_v, "velocity-y", filepath);
	}
}

void Cavern2DSimpleWorker::assemble_p_stroke_solver(){
	LodMatrix mat(p_size());
	for (size_t j = 0; j < _cc_grid.ny()+1; ++j)
	for (size_t i = 0; i < _cc_grid.nx()+1; ++i){
		bool is_left = (i==0);
		bool is_right = (i==_cc_grid.nx());
		bool is_bottom = (j==0);
		bool is_top = (j==_cc_grid.ny());

		size_t ind0 = _grid.cell_centered_grid_index_ip_jp(i, j);
		double coef_x = _du/_hx/_hx;
		double coef_y = _dv/_hy/_hy;
		// x
		if (!is_right){
			size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i+1, j);
			mat.add_value(ind0, ind0, coef_x);
			mat.add_value(ind0, ind1, -coef_x);
		}
		if (!is_left){
			size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i-1, j);
			mat.add_value(ind0, ind0, coef_x);
			mat.add_value(ind0, ind1, -coef_x);
		}
		// y
		if (!is_top){
			size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i, j+1);
			mat.add_value(ind0, ind0, coef_y);
			mat.add_value(ind0, ind1, -coef_y);
		}
		if (!is_bottom){
			size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i, j-1);
			mat.add_value(ind0, ind0, coef_y);
			mat.add_value(ind0, ind1, -coef_y);
		}
	}
	mat.set_unit_row(0);
	_p_stroke_solver.set_matrix(mat.to_csr());
}

void Cavern2DSimpleWorker::assemble_u_slae(){
	_rhs_u.resize(_u.size());
	std::fill(_rhs_u.begin(), _rhs_u.end(), 0.0);
	LodMatrix mat(_u.size());

	auto add_to_mat = [&](size_t row_index, std::array<size_t, 2> ij_col, double value){
		if (ij_col[1] == _grid.ny()){
			// ghost index => top boundary condition: u = 1
			size_t ind1 = _grid.yface_grid_index_i_jp(ij_col[0], ij_col[1]-1);
			mat.add_value(row_index, ind1, -value);
			_rhs_u[row_index] -= 2.0*value;
		} else if (ij_col[1] == (size_t)-1){
			// ghost index => bottom boundary condition: u = 0
			size_t ind1 = _grid.yface_grid_index_i_jp(ij_col[0], ij_col[1]+1);
			mat.add_value(row_index, ind1, -value);
			///////////////////////////////////////////////
			if (_u_bottom != 0) { _rhs_u[row_index] -= 2.0 * value; } // for sym cavern
			//////////////////////////////////////////////
		} else {
			size_t ind1 = _grid.yface_grid_index_i_jp(ij_col[0], ij_col[1]);
			mat.add_value(row_index, ind1, value);
		}
	};

	// left/right boundary: u = 0
	for (size_t j=0; j< _grid.ny(); ++j){
		size_t index_left = _grid.yface_grid_index_i_jp(0, j);
		add_to_mat(index_left, {0, j}, 1.0);
		_rhs_u[index_left] = 0.0;

		size_t index_right = _grid.yface_grid_index_i_jp(_grid.nx(), j);
		add_to_mat(index_right, {_grid.nx(), j}, 1.0);
		_rhs_u[index_right] = 0.0;
	}

	// internal
	for (size_t j=0; j < _grid.ny(); ++j)
	for (size_t i=1; i < _grid.nx(); ++i){
		size_t row_index = _grid.yface_grid_index_i_jp(i, j);   //[i, j+1/2] linear index in u grid

		double u0_plus   = u_ip_jp(i, j);   //_u[i+1/2, j+1/2]
		double u0_minus  = u_ip_jp(i-1, j); //_u[i-1/2, j+1/2]
		double v0_plus   = v_i_j(i, j+1);   // _v[i,j+1]
		double v0_minus  = v_i_j(i, j);     // _v[i,j]

		// u_(i,j+1/2)
		add_to_mat(row_index, {i, j}, 1.0);
		//     + tau * d(u0*u)/ dx
		add_to_mat(row_index, {i+1,j}, _tau/2.0/_hx*u0_plus);
		add_to_mat(row_index, {i-1,j}, -_tau/2.0/_hx*u0_minus);
		//     + tau * d(v0*u)/dy
		add_to_mat(row_index, {i, j+1}, _tau/2.0/_hy*v0_plus);
		add_to_mat(row_index, {i, j-1}, -_tau/2.0/_hy*v0_minus);
		//     - tau / Re * d^2u/dx^2
		add_to_mat(row_index, {i, j}, 2.0*_tau/_Re/_hx/_hy);
		add_to_mat(row_index, {i+1, j}, -_tau/_Re/_hx/_hx);
		add_to_mat(row_index, {i-1, j}, -_tau/_Re/_hy/_hy);
		//     - tau / Re * d^2u/dy^2
		add_to_mat(row_index, {i, j}, 2.0*_tau/_Re/_hy/_hy);
		add_to_mat(row_index, {i, j+1}, -_tau/_Re/_hy/_hy);
		add_to_mat(row_index, {i, j-1}, -_tau/_Re/_hy/_hy);
		// = u0_(i,j+1/2)
		_rhs_u[row_index] += _u[row_index];
		//      - tau * dp/dx
		_rhs_u[row_index] -= _tau/_hx*(p_ip_jp(i, j) - p_ip_jp(i-1, j));
	}
	_mat_u = mat.to_csr();
}

void Cavern2DSimpleWorker::assemble_v_slae(){
	_rhs_v.resize(_v.size());
	std::fill(_rhs_v.begin(), _rhs_v.end(), 0.0);
	LodMatrix mat(_v.size());

	auto add_to_mat = [&](size_t row_index, std::array<size_t, 2> ij_col, double value){
		if (ij_col[0] == (size_t)-1){
			// left boundary condition: v = 0
			size_t ind1 = _grid.xface_grid_index_ip_j(ij_col[0]+1, ij_col[1]);
			mat.add_value(row_index, ind1, -value);
		} else if (ij_col[0] == _grid.nx()){
			// right boundary condition: v = 0
			size_t ind1 = _grid.xface_grid_index_ip_j(ij_col[0]-1, ij_col[1]);
			mat.add_value(row_index, ind1, -value);
		} else {
			// internal
			size_t ind1 = _grid.xface_grid_index_ip_j(ij_col[0], ij_col[1]);
			mat.add_value(row_index, ind1, value);
		}
	};

	// === top/bottom boundaries
	for (size_t i=0; i<_grid.nx(); ++i){
		// top boundary condition: v = 0;
		size_t index_top = _grid.xface_grid_index_ip_j(i, _grid.ny());
		add_to_mat(index_top, {i, _grid.ny()}, 1.0);
		_rhs_v[index_top] = 0.0;
		// bottom boundary condition: v = 0;
		size_t index_bottom = _grid.xface_grid_index_ip_j(i, 0);
		add_to_mat(index_bottom, {i, 0}, 1.0);
		_rhs_v[index_bottom] = 0.0;
	}

	// === internal
	for (size_t j=1; j<_grid.ny(); ++j)
	for (size_t i=0; i<_grid.nx(); ++i){
		size_t row_index = _grid.xface_grid_index_ip_j(i, j);   //[i+1/2, j] linear index in v grid

		double u0_plus   = u_i_j(i+1, j);   // _u[i+1, j]
		double u0_minus  = u_i_j(i, j);     // _u[i, j]
		double v0_plus   = v_ip_jp(i, j);   // _v[i+1/2, j+1/2]
		double v0_minus  = v_ip_jp(i, j-1); // _v[i+1/2, j-1/2]

		// v_(i+1/2, j)
		add_to_mat(row_index, {i, j}, 1.0);
		//     + tau * d (u0*v) / dx
		add_to_mat(row_index, {i+1, j}, _tau/2.0/_hx*u0_plus);
		add_to_mat(row_index, {i-1, j}, -_tau/2.0/_hx*u0_minus);
		//     + tau * d (v0*v) / dy
		add_to_mat(row_index, {i, j+1}, _tau/2.0/_hy*v0_plus);
		add_to_mat(row_index, {i, j-1}, -_tau/2.0/_hy*v0_minus);
		//     - tau/Re * d^2v/dx^2
		add_to_mat(row_index, {i, j}, 2.0*_tau/_Re/_hx/_hx);
		add_to_mat(row_index, {i+1, j}, -_tau/_Re/_hx/_hx);
		add_to_mat(row_index, {i-1, j}, -_tau/_Re/_hx/_hx);
		//     - tau/Re * d^2v/dy^2
		add_to_mat(row_index, {i, j}, 2.0*_tau/_Re/_hy/_hy);
		add_to_mat(row_index, {i, j+1}, -_tau/_Re/_hy/_hy);
		add_to_mat(row_index, {i, j-1}, -_tau/_Re/_hy/_hy);
		// = v_prev
		_rhs_v[row_index] += _v[row_index];
		//     - tau * dp/dy
		_rhs_v[row_index] -= _tau*(p_ip_jp(i, j) - p_ip_jp(i, j-1))/_hy;
	}
	_mat_v = mat.to_csr();
}


std::vector<double> Cavern2DSimpleWorker::compute_u_star(){
	std::vector<double> u_star(_u);
	AmgcMatrixSolver::solve_slae(_mat_u, _rhs_u, u_star);
	return u_star;
}

std::vector<double> Cavern2DSimpleWorker::compute_v_star(){
	std::vector<double> v_star(_v);
	AmgcMatrixSolver::solve_slae(_mat_v, _rhs_v, v_star);
	return v_star;
}

std::vector<double> Cavern2DSimpleWorker::compute_p_stroke(const std::vector<double>& u_star, const std::vector<double>& v_star){
	std::vector<double> rhs(_p.size(), 0.0);
	for (size_t i = 0; i < _grid.nx(); ++i)
	for (size_t j = 0; j < _grid.ny(); ++j){
		size_t ind0 = _grid.cell_centered_grid_index_ip_jp(i, j);
		size_t ind_left = _grid.yface_grid_index_i_jp(i, j);
		size_t ind_right = _grid.yface_grid_index_i_jp(i+1, j);
		size_t ind_bot = _grid.xface_grid_index_ip_j(i, j);
		size_t ind_top = _grid.xface_grid_index_ip_j(i, j+1);
		rhs[ind0] = -(u_star[ind_right] - u_star[ind_left])/_tau/_hx - (v_star[ind_top] - v_star[ind_bot])/_tau/_hy;
	}
	rhs[0] = 0;
	std::vector<double> p_stroke;
	_p_stroke_solver.solve(rhs, p_stroke);
	return p_stroke;
}

std::vector<double> Cavern2DSimpleWorker::compute_u_stroke(const std::vector<double>& p_stroke){
	std::vector<double> u_stroke(_u.size(), 0.0);
	for (size_t i=1; i<_grid.nx(); ++i)
	for (size_t j=0; j<_grid.ny(); ++j){
		size_t ind0 = _grid.yface_grid_index_i_jp(i, j);
		size_t ind_plus  = _grid.cell_centered_grid_index_ip_jp(i, j);
		size_t ind_minus = _grid.cell_centered_grid_index_ip_jp(i-1, j);
		u_stroke[ind0] = -_tau * _du * (p_stroke[ind_plus] - p_stroke[ind_minus])/_hx;
	}
	return u_stroke;
}

std::vector<double> Cavern2DSimpleWorker::compute_v_stroke(const std::vector<double>& p_stroke){
	std::vector<double> v_stroke(_v.size(), 0.0);
	for (size_t i=0; i<_grid.nx(); ++i)
	for (size_t j=1; j<_grid.ny(); ++j){
		size_t ind0 = _grid.xface_grid_index_ip_j(i, j);
		size_t ind_plus  = _grid.cell_centered_grid_index_ip_jp(i, j);
		size_t ind_minus = _grid.cell_centered_grid_index_ip_jp(i, j-1);
		v_stroke[ind0] = -_tau * _dv * (p_stroke[ind_plus] - p_stroke[ind_minus])/_hy;
	}
	return v_stroke;
}

double Cavern2DSimpleWorker::u_i_j(size_t i, size_t j) const{
	size_t ind0 = _grid.yface_grid_index_i_jp(i, j);
	size_t ind1 = _grid.yface_grid_index_i_jp(i, j-1);
	return (_u[ind0] + _u[ind1])/2.0;
}
double Cavern2DSimpleWorker::u_ip_jp(size_t i, size_t j) const{
	size_t ind0 = _grid.yface_grid_index_i_jp(i, j);
	size_t ind1 = _grid.yface_grid_index_i_jp(i+1, j);
	return (_u[ind0] + _u[ind1])/2.0;
}

double Cavern2DSimpleWorker::v_i_j(size_t i, size_t j) const{
	size_t ind0 = _grid.xface_grid_index_ip_j(i, j);
	size_t ind1 = _grid.xface_grid_index_ip_j(i-1, j);
	return (_v[ind0] + _v[ind1])/2.0;
}
double Cavern2DSimpleWorker::v_ip_jp(size_t i, size_t j) const{
	size_t ind0 = _grid.xface_grid_index_ip_j(i, j);
	size_t ind1 = _grid.xface_grid_index_ip_j(i, j+1);
	return (_v[ind0] + _v[ind1])/2.0;
}

double Cavern2DSimpleWorker::p_ip_jp(size_t i, size_t j) const{
	return _p[_grid.cell_centered_grid_index_ip_jp(i, j)];
}

std::vector<Vector> Cavern2DSimpleWorker::build_main_grid_velocity() const{
	std::vector<Vector> ret(_grid.n_points());
	// boundary
	for (size_t j = 0; j < _grid.ny()+1; ++j){
		// left boundary
		size_t ind_left = _grid.to_linear_point_index({0, j});
		ret[ind_left] = Vector(0, 0);
		// right boundary
		size_t ind_right = _grid.to_linear_point_index({_grid.nx(), j});
		ret[ind_right] = Vector(0, 0);
	}
	for (size_t i = 0; i < _grid.nx()+1; ++i){
		// bottom boundary
		size_t ind_bot = _grid.to_linear_point_index({i, 0});
		ret[ind_bot] = Vector(get_bottom_velocity(), 0); ////////////////////////////////////// ��� ��� 0 //////////////////////////
		// top boundary
		size_t ind_top = _grid.to_linear_point_index({i, _grid.ny()});
		ret[ind_top] = Vector(1, 0);
	}
	// internal
	for (size_t j=1; j<_grid.ny(); ++j)
	for (size_t i=1; i<_grid.nx(); ++i){
		size_t ind = _grid.to_linear_point_index({i, j});
		size_t ind_top = _grid.yface_grid_index_i_jp(i, j);
		size_t ind_bot = _grid.yface_grid_index_i_jp(i, j-1);
		size_t ind_left = _grid.xface_grid_index_ip_j(i-1, j);
		size_t ind_right = _grid.xface_grid_index_ip_j(i, j);
		ret[ind] = Vector(0.5*(_u[ind_top] + _u[ind_bot]),
				  0.5*(_v[ind_left] + _v[ind_right]));
	}

	return ret;
}


//////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Cavern 2D, SIMPLE algorithm", "[cavern2-simple]"){
	std::cout << std::endl << "--- cfd24_test [cavern2-simple] --- " << std::endl;

	norms.open("C:/git_repos/CFDCourse24/build/src/test/norms.txt");

	// problem parameters
	double Re = 100;
	double tau = 0.01;
	double alpha = 0.08;
	size_t n_cells = 30;
	size_t max_it = 1000;
	double eps = 0.01;
	double u_bottom = 0.0;

	// worker initialization
	Cavern2DSimpleWorker worker(Re, n_cells, tau, alpha, u_bottom);
	worker.initialize_saver(false, "cavern2");

	// initial condition
	std::vector<double> u_init(worker.u_size(), 0.0);
	std::vector<double> v_init(worker.v_size(), 0.0);
	std::vector<double> p_init(worker.p_size(), 0.0);
	worker.set_uvp(u_init, v_init, p_init);
	worker.save_current_fields(0);

	// iterations loop
	size_t it = 0;
	for (it=1; it < max_it; ++it){
		double nrm = worker.step();

		// print norm and pressure value at the top-right corner
		std::cout << it << " " << nrm << " " << worker.pressure().back() << std::endl;

		// export solution to vtk
		worker.save_current_fields(it);

		// break if residual is low enough
		if (nrm < eps){
			break;
		}
	}
	CHECK(it == 9);

	norms.close();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Cavern2DSimpleSymWorker : Cavern2DSimpleWorker {
	Cavern2DSimpleSymWorker(double Re, size_t n_cells, double tau, double alpha_p, double u_bottom) : 
		Cavern2DSimpleWorker(Re, n_cells, tau, alpha_p, u_bottom) {}

private:
	//void assemble_u_slae();
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Cavern 2D, SIMPLE algorithm sym", "[cavern2-simple-sym]") {
	std::cout << std::endl << "--- cfd24_test [cavern2-simple-sym] --- " << std::endl;

	// problem parameters
	double Re = 100;
	double tau = 0.03;
	double alpha = 0.8;
	size_t n_cells = 30;
	size_t max_it = 1000;
	double eps = 1e-1;
	double u_bottom = 1.0;

	// worker initialization
	Cavern2DSimpleSymWorker worker(Re, n_cells, tau, alpha, u_bottom);
	worker.initialize_saver(false, "cavern2sym");

	// initial condition
	std::vector<double> u_init(worker.u_size(), 0.0);
	std::vector<double> v_init(worker.v_size(), 0.0);
	std::vector<double> p_init(worker.p_size(), 0.0);
	worker.set_uvp(u_init, v_init, p_init);
	worker.save_current_fields(0);

	// iterations loop
	size_t it = 0;
	for (it = 1; it < max_it; ++it) {
		double nrm = worker.step();

		// print norm and pressure value at the top-right corner
		std::cout << it << " " << nrm << " " << worker.pressure().back() << std::endl;

		// export solution to vtk
		worker.save_current_fields(it);

		// break if residual is low enough
		if (nrm < eps) {
			break;
		}
	}
	CHECK(it == 9);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////