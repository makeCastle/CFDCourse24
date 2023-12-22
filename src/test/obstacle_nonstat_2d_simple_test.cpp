#include "cfd24_test.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/mat/lodmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/debug/printer.hpp"
#include "utils/vecmat.hpp"
#include "cfd24/debug/tictoc.hpp"
#include <sstream>
#include <iomanip>

using namespace cfd;

struct ObstacleNonstat2DSimpleWorker{
	ObstacleNonstat2DSimpleWorker(double Re, double Pe, const RegularGrid2D& grid, double E, double time_step);

	void initialize();
	double set_uvp(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& p);
	double step();

	void initialize_saver(bool save_exact_fields, std::string stem);
	void save_current_fields(double time);

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
	double to_next_time_step();

	struct Coefficients{
		double cpx1, cpx2;
		double cpy1, cpy2;
		double cfx1, cfx2;
		double cfy1, cfy2;
		double cx1, cx2;
		double cy1, cy2;
		double nu1, nu2;
	};
	Coefficients coefficients() const;


private:
	const RegularGrid2D _grid;
	const RegularGrid2D _cc_grid;
	const RegularGrid2D _xf_grid;
	const RegularGrid2D _yf_grid;
	const double _hx;
	const double _hy;
	const double _Re;
	const double _Pe;
	const double _tau;
	const double _alpha_p;
	const double _time_step;

	std::vector<double> _p;
	std::vector<double> _u;
	std::vector<double> _v;
	std::vector<double> _t;
	std::vector<double> _u_old;
	std::vector<double> _v_old;

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
	std::ofstream cx_writer;

	void assemble_p_stroke_solver();
	void assemble_u_slae();
	void assemble_v_slae();

	std::vector<double> compute_phi() const;
	std::vector<double> compute_u_stroke_outflow(const std::vector<double>& u_star) const;
	std::vector<double> compute_u_star();
	std::vector<double> compute_v_star();
	std::vector<double> compute_p_stroke(const std::vector<double>& u_star, const std::vector<double>& v_star, const std::vector<double>& u_stroke_outflow);
	std::vector<double> compute_u_stroke(const std::vector<double>& p_stroke, const std::vector<double>& u_stroke_outflow);
	std::vector<double> compute_v_stroke(const std::vector<double>& p_stroke);
	std::vector<double> compute_temperature() const;
	double u_i_j(size_t i, size_t j) const;
	double u_ip_jp(size_t i, size_t j) const;
	double v_i_j(size_t i, size_t j) const;
	double v_ip_jp(size_t i, size_t j) const;
	double p_ip_jp(size_t i, size_t j) const;
	std::vector<Vector> build_main_grid_velocity() const;



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int gamma_closest_to_cell(size_t icell) const;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

};

ObstacleNonstat2DSimpleWorker::ObstacleNonstat2DSimpleWorker(double Re, double Pe, const RegularGrid2D& grid, double E, double time_step):
	_grid(grid),
	_cc_grid(_grid.cell_centered_grid()),
	_xf_grid(_grid.xface_centered_grid()),
	_yf_grid(_grid.yface_centered_grid()),
	_hx(_grid.Lx()/_grid.nx()),
	_hy(_grid.Ly()/_grid.ny()),
	_Re(Re),
	_Pe(Pe),
	_tau(E*Re*harmonic_mean(_hx*_hx, _hy*_hy)/4.0),
	_alpha_p(1.0/(1.0 + E)),
	_time_step(time_step)
{
	_du = 1.0 / (1 + _tau/_time_step + 2.0*_tau/_Re * (1.0/_hx/_hx + 1.0/_hy/_hy));
	_dv = 1.0 / (1 + _tau/_time_step + 2.0*_tau/_Re * (1.0/_hx/_hx + 1.0/_hy/_hy));
	assemble_p_stroke_solver();
}

std::vector<double> ObstacleNonstat2DSimpleWorker::compute_phi() const{
	cfd::LodMatrix lhs(_grid.n_cells());
	std::vector<double> rhs(_grid.n_cells(), 0);
	double coef_x = _hy/_hx;
	double coef_y = _hx/_hy;
	for (size_t j = 0; j < _grid.ny(); ++j)
	for (size_t i = 0; i < _grid.nx(); ++i){
		size_t ind_center = _grid.cell_centered_grid_index_ip_jp(i, j);
		if (_grid.is_active_cell(ind_center)){
			if (i > 0){
				size_t ind_left = _grid.cell_centered_grid_index_ip_jp(i-1, j);
				if (_grid.is_active_cell(ind_left)){
					lhs.add_value(ind_center, ind_center, coef_x);
					lhs.add_value(ind_center, ind_left, -coef_x);
				}
			} else {
				rhs[ind_center] -= _hy;
			}
			if (i < _grid.nx()-1){
				size_t ind_right = _grid.cell_centered_grid_index_ip_jp(i+1, j);
				if (_grid.is_active_cell(ind_right)){
					lhs.add_value(ind_center, ind_center, coef_x);
					lhs.add_value(ind_center, ind_right, -coef_x);
				}
			} else {
				rhs[ind_center] += _hy;
			}
			if (j > 0){
				size_t ind_bot = _grid.cell_centered_grid_index_ip_jp(i, j-1);
				if (_grid.is_active_cell(ind_bot)){
					lhs.add_value(ind_center, ind_center, coef_y);
					lhs.add_value(ind_center, ind_bot, -coef_y);
				}
			}
			if (j < _grid.ny()-1){
				size_t ind_top = _grid.cell_centered_grid_index_ip_jp(i, j+1);
				if (_grid.is_active_cell(ind_top)){
					lhs.add_value(ind_center, ind_center, coef_y);
					lhs.add_value(ind_center, ind_top, -coef_y);
				}
			}
		} else {
			lhs.add_value(ind_center, ind_center, 1.0);
		}
	}
	lhs.set_unit_row(0);
	rhs[0] = 0;
	std::vector<double> phi;
	AmgcMatrixSolver::solve_slae(lhs.to_csr(), rhs, phi);

	_cc_grid.save_vtk("phi.vtk");
	VtkUtils::add_point_data(phi, "phi", "phi.vtk");

	return phi;
}

void ObstacleNonstat2DSimpleWorker::initialize(){
	std::vector<double> phi = compute_phi();
	// build u
	std::vector<double> u(u_size(), 0);
	for (size_t j=0; j < _grid.ny(); ++j)
	for (size_t i=0; i < _grid.nx()+1; ++i){
		size_t ind0 = _grid.yface_grid_index_i_jp(i, j);
		if (i == 0 || i == _grid.nx()){
			u[ind0] = 1;
		} else {
			size_t ind_right = _grid.cell_centered_grid_index_ip_jp(i, j);
			size_t ind_left = _grid.cell_centered_grid_index_ip_jp(i-1, j);
			if (_grid.is_active_cell(ind_right) && _grid.is_active_cell(ind_left)){
				u[ind0] = (phi[ind_right] - phi[ind_left])/_hx;
			} else {
				u[ind0] = 0;
			}
		}
	}

	// build v
	std::vector<double> v(v_size(), 0);
	for (size_t j=0; j < _grid.ny()+1; ++j)
	for (size_t i=0; i < _grid.nx(); ++i){
		size_t ind0 = _grid.xface_grid_index_ip_j(i, j);
		if (j == 0 || j == _grid.ny()){
			v[ind0] = 0;
		} else {
			size_t ind_top = _grid.cell_centered_grid_index_ip_jp(i, j);
			size_t ind_bot = _grid.cell_centered_grid_index_ip_jp(i, j-1);
			if (_grid.is_active_cell(ind_top) && _grid.is_active_cell(ind_bot)){
				v[ind0] = (phi[ind_top] - phi[ind_bot])/_hy;
			} else {
				v[ind0] = 0;
			}
		}
	}

	_u = u;
	_v = v;
	_p = std::vector<double>(p_size(), 0);
	_t = std::vector<double>(p_size(), 0);

	to_next_time_step();
};

void ObstacleNonstat2DSimpleWorker::initialize_saver(bool save_exact_fields, std::string stem){
	_writer_all.reset(new VtkUtils::TimeDependentWriter(stem));
	if (save_exact_fields){
		_writer_u.reset(new VtkUtils::TimeDependentWriter(stem + "-u"));
		_writer_v.reset(new VtkUtils::TimeDependentWriter(stem + "-v"));
		_writer_p.reset(new VtkUtils::TimeDependentWriter(stem + "-p"));
	}
	cx_writer.open("c.txt");
	cx_writer << "Time Cx1 Cy1 Nu1 Cx2 Cy2 Nu2" << std::endl;

};

double ObstacleNonstat2DSimpleWorker::set_uvp(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& p){
	_u = u;
	_v = v;
	_p = p;
	assemble_u_slae();
	assemble_v_slae();
	// residuals
	double nrm_u = compute_residual(_mat_u, _rhs_u, _u)/_tau;
	double nrm_v = compute_residual(_mat_v, _rhs_v, _v)/_tau;
	return std::max(nrm_u, nrm_v);
};

double ObstacleNonstat2DSimpleWorker::to_next_time_step(){
	_t = compute_temperature();
	_u_old = _u;
	_v_old = _v;
	return set_uvp(_u, _v, _p);
}

double ObstacleNonstat2DSimpleWorker::step(){
	// Predictor step: U-star
	std::vector<double> u_star = compute_u_star();
	std::vector<double> v_star = compute_v_star();
	std::vector<double> u_stroke_outflow = compute_u_stroke_outflow(u_star);
	// Pressure correction
	std::vector<double> p_stroke = compute_p_stroke(u_star, v_star, u_stroke_outflow);
	// Velocity correction
	std::vector<double> u_stroke = compute_u_stroke(p_stroke, u_stroke_outflow);
	std::vector<double> v_stroke = compute_v_stroke(p_stroke);
	// Set final values
	std::vector<double> u_new = vector_sum(u_star, 1.0, u_stroke);
	std::vector<double> v_new = vector_sum(v_star, 1.0, v_stroke);
	std::vector<double> p_new = vector_sum(_p, _alpha_p, p_stroke);

	return set_uvp(u_new, v_new, p_new);
}

void ObstacleNonstat2DSimpleWorker::save_current_fields(double time){
	std::vector<double> actnum(_grid.actnum().begin(), _grid.actnum().end());
	// all data on the main grid
	if (_writer_all){
		std::string filepath = _writer_all->add(time);
		_grid.save_vtk(filepath);
		VtkUtils::add_cell_data(_p, "pressure", filepath);
		VtkUtils::add_cell_data(_t, "temperature", filepath);
		VtkUtils::add_cell_data(actnum, "actnum", filepath);
		VtkUtils::add_point_vector(build_main_grid_velocity(), "velocity", filepath);
	}
	// pressure
	if (_writer_p){
		std::string filepath = _writer_p->add(time);
		_cc_grid.save_vtk(filepath);
		VtkUtils::add_point_data(_p, "pressure", filepath);
		VtkUtils::add_point_data(_t, "temparature", filepath);
		VtkUtils::add_point_data(actnum, "actnum", filepath);
	}
	// u
	if (_writer_u){
		std::string filepath = _writer_u->add(time);
		_yf_grid.save_vtk(filepath);
		VtkUtils::add_point_data(_u, "velocity-x", filepath);
	}
	// v
	if (_writer_v){
		std::string filepath = _writer_v->add(time);
		_xf_grid.save_vtk(filepath);
		VtkUtils::add_point_data(_v, "velocity-y", filepath);
	}
	
	// write coefficients to file
	Coefficients coefs = coefficients();
	cx_writer << time << " ";
	cx_writer << coefs.cx1 << " " << coefs.cy1 << " " << coefs.nu1 << " ";
	cx_writer << coefs.cx2 << " " << coefs.cy2 << " " << coefs.nu2 << std::endl;
}


void ObstacleNonstat2DSimpleWorker::assemble_p_stroke_solver(){
	LodMatrix mat(p_size());
	double coef_x = _du/_hx/_hx;
	double coef_y = _dv/_hy/_hy;
	for (size_t j = 0; j < _cc_grid.ny()+1; ++j)
	for (size_t i = 0; i < _cc_grid.nx()+1; ++i){
		size_t ind0 = _grid.cell_centered_grid_index_ip_jp(i, j);
		if (_grid.is_active_cell(ind0)){
			bool is_left = (i==0);
			bool is_right = (i==_cc_grid.nx());
			bool is_bottom = (j==0);
			bool is_top = (j==_cc_grid.ny());
			// x
			if (!is_right){
				size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i+1, j);
				if (_grid.is_active_cell(ind1)){
					mat.add_value(ind0, ind0, coef_x);
					mat.add_value(ind0, ind1, -coef_x);
				}
			}
			if (!is_left){
				size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i-1, j);
				if (_grid.is_active_cell(ind1)){
					mat.add_value(ind0, ind0, coef_x);
					mat.add_value(ind0, ind1, -coef_x);
				}
			}
			// y
			if (!is_top){
				size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i, j+1);
				if (_grid.is_active_cell(ind1)){
					mat.add_value(ind0, ind0, coef_y);
					mat.add_value(ind0, ind1, -coef_y);
				}
			}
			if (!is_bottom){
				size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i, j-1);
				if (_grid.is_active_cell(ind1)){
					mat.add_value(ind0, ind0, coef_y);
					mat.add_value(ind0, ind1, -coef_y);
				}
			}
		} else {
			mat.set_value(ind0, ind0, 1.0);
		}
	}
	//regularization
	//double diagval = 1e-6*(coef_x + coef_y)/2;
	//for (size_t i=0; i<_cc_grid.n_points(); ++i){
		//mat.add_value(i, i, diagval);
	//}
	mat.set_unit_row(0);
	_p_stroke_solver.set_matrix(mat.to_csr());
}

void ObstacleNonstat2DSimpleWorker::assemble_u_slae(){
	_rhs_u.resize(_u.size());
	std::fill(_rhs_u.begin(), _rhs_u.end(), 0.0);
	LodMatrix mat(_u.size());

	auto add_to_mat = [&](size_t row_index, std::array<size_t, 2> ij_col, double value){
		if (ij_col[1] == _grid.ny()){
			// ghost index => top boundary condition: du/dn = 0
			size_t ind1 = _grid.yface_grid_index_i_jp(ij_col[0], ij_col[1]-1);
			mat.add_value(row_index, ind1, value);
		} else if (ij_col[1] == (size_t)-1){
			// ghost index => bottom boundary condition: du/dn = 0
			size_t ind1 = _grid.yface_grid_index_i_jp(ij_col[0], ij_col[1]+1);
			mat.add_value(row_index, ind1, value);
		} else {
			size_t ind1 = _grid.yface_grid_index_i_jp(ij_col[0], ij_col[1]);
			if (_grid.yface_type(ind1) == RegularGrid2D::FaceType::Deactivated){
				// ghost index => obstacle boundary u = 0
				mat.add_value(row_index, row_index, -value);
			} else {
				mat.add_value(row_index, ind1, value);
			}
		}
	};

	for (size_t j=0; j< _grid.ny(); ++j){
		// left boundary: u = 1
		{
			size_t index_left = _grid.yface_grid_index_i_jp(0, j);
			add_to_mat(index_left, {0, j}, 1.0);
			_rhs_u[index_left] = 1.0;
		}
		// right boundary: du/dt + u*du/dx = 0
		{
			size_t index_right = _grid.yface_grid_index_i_jp(_grid.nx(), j);
			size_t index_right_prev = _grid.yface_grid_index_i_jp(_grid.nx()-1, j);
			double u0 = std::max(0.0, _u[index_right]);
			double coef = _tau*u0/_hx;
			mat.set_value(index_right, index_right, 1.0 + _tau/_time_step + coef);
			mat.set_value(index_right, index_right_prev, -coef);
			_rhs_u[index_right] = u0 + _tau/_time_step * _u_old[index_right];
		}
	}

	// internal
	for (size_t j=0; j < _grid.ny(); ++j)
	for (size_t i=1; i < _grid.nx(); ++i){
		size_t row_index = _grid.yface_grid_index_i_jp(i, j);   //[i, j+1/2] linear index in u grid
		if (_grid.yface_type(row_index) == RegularGrid2D::FaceType::Internal){
			double u0_plus   = u_ip_jp(i, j);   // _u[i+1/2, j+1/2]
			double u0_minus  = u_ip_jp(i-1, j); // _u[i-1/2, j+1/2]
			double v0_plus   = v_i_j(i, j+1);   // _v[i,j+1]
			double v0_minus  = v_i_j(i, j);     // _v[i,j]

			// (1 + tau/time_step) * u_(i,j+1/2)
			add_to_mat(row_index, {i, j}, 1.0 + _tau/_time_step);
			//     + tau * d(u0*u)/ dx
			add_to_mat(row_index, {i+1,j}, _tau/2.0/_hx*u0_plus);
			add_to_mat(row_index, {i-1,j}, -_tau/2.0/_hx*u0_minus);
			//     + tau * d(v0*u)/dy
			add_to_mat(row_index, {i, j+1}, _tau/2.0/_hy*v0_plus);
			add_to_mat(row_index, {i, j-1}, -_tau/2.0/_hy*v0_minus);
			//     - tau / Re * d^2u/dx^2
			add_to_mat(row_index, {i, j}, 2.0*_tau/_Re/_hx/_hx);
			add_to_mat(row_index, {i+1, j}, -_tau/_Re/_hx/_hx);
			add_to_mat(row_index, {i-1, j}, -_tau/_Re/_hx/_hx);
			//     - tau / Re * d^2u/dy^2
			add_to_mat(row_index, {i, j}, 2.0*_tau/_Re/_hy/_hy);
			add_to_mat(row_index, {i, j+1}, -_tau/_Re/_hy/_hy);
			add_to_mat(row_index, {i, j-1}, -_tau/_Re/_hy/_hy);
			// = u0_(i,j+1/2)
			_rhs_u[row_index] += _u[row_index];
			// + tau/time_step*u_old(i, j+1/2)
			_rhs_u[row_index] += (_tau/_time_step)*_u_old[row_index];
			//      - tau * dp/dx
			_rhs_u[row_index] -= _tau/_hx*(p_ip_jp(i, j) - p_ip_jp(i-1, j));
		} else {
			mat.set_value(row_index, row_index, 1.0);
			_rhs_u[row_index] = 0;
		}
	}
	_mat_u = mat.to_csr();
}

void ObstacleNonstat2DSimpleWorker::assemble_v_slae(){
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
			if (_grid.xface_type(ind1) == RegularGrid2D::FaceType::Deactivated){
				mat.add_value(row_index, row_index, -value);
			} else {
				mat.add_value(row_index, ind1, value);
			}
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
	
		if (_grid.xface_type(row_index) == RegularGrid2D::FaceType::Internal){
			double u0_plus   = u_i_j(i+1, j);   // _u[i+1, j]
			double u0_minus  = u_i_j(i, j);     // _u[i, j]
			double v0_plus   = v_ip_jp(i, j);   // _v[i+1/2, j+1/2]
			double v0_minus  = v_ip_jp(i, j-1); // _v[i+1/2, j-1/2]

			// (1 + gamma/tau) * v_(i+1/2, j)
			add_to_mat(row_index, {i, j}, 1.0 + _tau/_time_step);
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
			// = v
			_rhs_v[row_index] += _v[row_index];
			//   + (gamma/tau)*v_old
			_rhs_v[row_index] += (_tau/_time_step)*_v_old[row_index];
			//     - tau * dp/dy
			_rhs_v[row_index] -= _tau*(p_ip_jp(i, j) - p_ip_jp(i, j-1))/_hy;
		} else {
			// v = 0 on boundary and deactivated faces
			mat.set_value(row_index, row_index, 1.0);
			_rhs_v[row_index] = 0;
		}
	}
	_mat_v = mat.to_csr();
}


std::vector<double> ObstacleNonstat2DSimpleWorker::compute_u_star(){
	std::vector<double> u_star(_u);
	AmgcMatrixSolver::solve_slae(_mat_u, _rhs_u, u_star);
	return u_star;
}

std::vector<double> ObstacleNonstat2DSimpleWorker::compute_v_star(){
	std::vector<double> v_star(_v);
	AmgcMatrixSolver::solve_slae(_mat_v, _rhs_v, v_star);
	return v_star;
}

std::vector<double> ObstacleNonstat2DSimpleWorker::compute_p_stroke(
		const std::vector<double>& u_star,
		const std::vector<double>& v_star,
		const std::vector<double>& u_stroke_outflow){
	std::vector<double> rhs(_p.size(), 0.0);
	for (size_t i = 0; i < _grid.nx(); ++i)
	for (size_t j = 0; j < _grid.ny(); ++j){
		size_t ind0 = _grid.cell_centered_grid_index_ip_jp(i, j);
		size_t ind_left = _grid.yface_grid_index_i_jp(i, j);
		size_t ind_right = _grid.yface_grid_index_i_jp(i+1, j);
		size_t ind_bot = _grid.xface_grid_index_ip_j(i, j);
		size_t ind_top = _grid.xface_grid_index_ip_j(i, j+1);
		rhs[ind0] = -(u_star[ind_right] - u_star[ind_left])/_tau/_hx - (v_star[ind_top] - v_star[ind_bot])/_tau/_hy;
		// outflow compensation
		if (i == _grid.nx()-1){
			rhs[ind0] -= (u_stroke_outflow[j]) / _tau / _hx;
		}
	}
	rhs[0] = 0;
	std::vector<double> p_stroke;
	_p_stroke_solver.solve(rhs, p_stroke);
	double p0 = (p_stroke[0] + p_stroke[_grid.cell_centered_grid_index_ip_jp(0, _grid.ny()-1)])/2.0;
	for (size_t i=0; i<p_stroke.size(); ++i){
		if (_grid.is_active_cell(i)){
			p_stroke[i] -= p0;
		} else {
			p_stroke[i] = 0;
		}
	}
	return p_stroke;
}

std::vector<double> ObstacleNonstat2DSimpleWorker::compute_u_stroke_outflow(const std::vector<double>& u_star) const{
	double qin = 0;
	double qout = 0;
	for (size_t j=0; j<_grid.ny(); ++j){
		size_t ind_left = _grid.yface_grid_index_i_jp(0, j);
		size_t ind_right = _grid.yface_grid_index_i_jp(_grid.nx(), j);
		qin += u_star[ind_left]*_hy;
		qout += u_star[ind_right]*_hy;
	}
	double Lout = _grid.Ly();
	double diff_u = (qin - qout)/Lout;
	std::vector<double> ret(_grid.ny(), diff_u);
	return ret;
}

std::vector<double> ObstacleNonstat2DSimpleWorker::compute_u_stroke(const std::vector<double>& p_stroke, const std::vector<double>& u_stroke_outflow){
	std::vector<double> u_stroke(_u.size(), 0.0);
	for (size_t i=1; i<_grid.nx(); ++i)
	for (size_t j=0; j<_grid.ny(); ++j){
		size_t ind0 = _grid.yface_grid_index_i_jp(i, j);
		if (_grid.yface_type(ind0) == RegularGrid2D::FaceType::Internal){
			size_t ind_plus  = _grid.cell_centered_grid_index_ip_jp(i, j);
			size_t ind_minus = _grid.cell_centered_grid_index_ip_jp(i-1, j);
			u_stroke[ind0] = -_tau * _du * (p_stroke[ind_plus] - p_stroke[ind_minus])/_hx;
		}
	}
	// outflow
	for (size_t j=0; j<_grid.ny(); ++j){
		size_t ind0 = _grid.yface_grid_index_i_jp(_grid.nx(), j);
		u_stroke[ind0] = u_stroke_outflow[j];
	}
	return u_stroke;
}

std::vector<double> ObstacleNonstat2DSimpleWorker::compute_v_stroke(const std::vector<double>& p_stroke){
	std::vector<double> v_stroke(_v.size(), 0.0);
	for (size_t i=0; i<_grid.nx(); ++i)
	for (size_t j=1; j<_grid.ny(); ++j){
		size_t ind0 = _grid.xface_grid_index_ip_j(i, j);
		if (_grid.xface_type(ind0) == RegularGrid2D::FaceType::Internal){
			size_t ind_plus  = _grid.cell_centered_grid_index_ip_jp(i, j);
			size_t ind_minus = _grid.cell_centered_grid_index_ip_jp(i, j-1);
			v_stroke[ind0] = -_tau * _dv * (p_stroke[ind_plus] - p_stroke[ind_minus])/_hy;
		}
	}
	return v_stroke;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ObstacleNonstat2DSimpleWorker::gamma_closest_to_cell(size_t icell) const {
	double x = _grid.cell_center(icell).x();
	if (x < 3.25) { // center between first and second obstacles
		return 1;
	}
	else {
		return 2;
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double> ObstacleNonstat2DSimpleWorker::compute_temperature() const{
	std::vector<double> rhs(_p.size(), 0);
	LodMatrix mat(_p.size());
	auto add_to_mat = [&](size_t row_index, std::array<size_t, 2> ij_col, double value){
		if (ij_col[0] == (size_t)-1){
			// left boundary: T=0
			mat.add_value(row_index, row_index, -value);
		} else if (ij_col[0] == _grid.nx()){
			// right boundary: dT/dn = 0
			mat.add_value(row_index, row_index, value);
		} else if (ij_col[1] == (size_t)-1){
			// bottom boundary: dT/dn = 0
			mat.add_value(row_index, row_index, value);
		} else if (ij_col[1] == _grid.ny()){
			// bottom boundary: dT/dn = 0
			mat.add_value(row_index, row_index, value);
		} else {
			size_t ind1 = _grid.cell_centered_grid_index_ip_jp(ij_col[0], ij_col[1]);
			if (_grid.is_active_cell(ind1)){
				mat.add_value(row_index, ind1, value);
			} else {
				// internal boundary: T = 1
				mat.add_value(row_index, row_index, -value);
				double t_gamma = (gamma_closest_to_cell(row_index) == 1) ? 0.5 : 1.0;
				rhs[row_index] -= 2 * t_gamma * value;
			}
		}
	};
	for (size_t j=0; j < _grid.ny(); ++j)
	for (size_t i=0; i < _grid.nx(); ++i){
		size_t row_index = _grid.cell_centered_grid_index_ip_jp(i, j);
		if (_grid.is_active_cell(row_index)){
			double u_left = _u[_grid.yface_grid_index_i_jp(i, j)];
			double u_right = _u[_grid.yface_grid_index_i_jp(i+1, j)];
			double v_bot = _v[_grid.xface_grid_index_ip_j(i, j)];
			double v_top = _v[_grid.xface_grid_index_ip_j(i, j+1)];

			// 1.0/time_step T(i+1/2,j+1/2)
			add_to_mat(row_index, {i, j}, 1.0 / _time_step);
			//     + d(u0*T)/ dx
			add_to_mat(row_index, {i+1,j},  u_right/2.0/_hx);
			add_to_mat(row_index, {i-1,j}, -u_left/2.0/_hx);
			//     + d(v0*T)/dy
			add_to_mat(row_index, {i, j+1},  v_top/2.0/_hy);
			add_to_mat(row_index, {i, j-1}, -v_bot/2.0/_hy);
			//     - 1 / Re * d^2u/dx^2
			add_to_mat(row_index, {i, j},    2.0/_Pe/_hx/_hx);
			add_to_mat(row_index, {i+1, j}, -1.0/_Pe/_hx/_hx);
			add_to_mat(row_index, {i-1, j}, -1.0/_Pe/_hx/_hx);
			//     - 1 / Re * d^2u/dy^2
			add_to_mat(row_index, {i, j},     2.0/_Pe/_hy/_hy);
			add_to_mat(row_index, {i, j+1},  -1.0/_Pe/_hy/_hy);
			add_to_mat(row_index, {i, j-1},  -1.0/_Pe/_hy/_hy);
			// = 1.0 / time_step*Told
			rhs[row_index] += 1.0 / _time_step*_t[row_index];
		} else {
			mat.set_value(row_index, row_index, 1.0);
			rhs[row_index] = 0;
		}
	}
	std::vector<double> temperature;
	AmgcMatrixSolver::solve_slae(mat.to_csr(), rhs, temperature);
	return temperature;
}

double ObstacleNonstat2DSimpleWorker::u_i_j(size_t i, size_t j) const{
	size_t ind0 = _grid.yface_grid_index_i_jp(i, j);
	size_t ind1 = _grid.yface_grid_index_i_jp(i, j-1);
	return (_u[ind0] + _u[ind1])/2.0;
}
double ObstacleNonstat2DSimpleWorker::u_ip_jp(size_t i, size_t j) const{
	size_t ind0 = _grid.yface_grid_index_i_jp(i, j);
	size_t ind1 = _grid.yface_grid_index_i_jp(i+1, j);
	return (_u[ind0] + _u[ind1])/2.0;
}

double ObstacleNonstat2DSimpleWorker::v_i_j(size_t i, size_t j) const{
	size_t ind0 = _grid.xface_grid_index_ip_j(i, j);
	size_t ind1 = _grid.xface_grid_index_ip_j(i-1, j);
	return (_v[ind0] + _v[ind1])/2.0;
}
double ObstacleNonstat2DSimpleWorker::v_ip_jp(size_t i, size_t j) const{
	size_t ind0 = _grid.xface_grid_index_ip_j(i, j);
	size_t ind1 = _grid.xface_grid_index_ip_j(i, j+1);
	return (_v[ind0] + _v[ind1])/2.0;
}

double ObstacleNonstat2DSimpleWorker::p_ip_jp(size_t i, size_t j) const{
	return _p[_grid.cell_centered_grid_index_ip_jp(i, j)];
}

std::vector<Vector> ObstacleNonstat2DSimpleWorker::build_main_grid_velocity() const{
	std::vector<Vector> ret(_grid.n_points());
	// boundary
	for (size_t i = 0; i < _grid.nx()+1; ++i){
		{
			// bottom boundary
			size_t ind_bot = _grid.to_linear_point_index({i, 0});
			double u = _u[_grid.yface_grid_index_i_jp(i, 0)];
			ret[ind_bot] = Vector(u, 0);
		}
		{
			// top boundary
			size_t ind_top = _grid.to_linear_point_index({i, _grid.ny()});
			double u = _u[_grid.yface_grid_index_i_jp(i, _grid.ny()-1)];
			ret[ind_top] = Vector(u, 0);
		}
	}
	for (size_t j = 0; j < _grid.ny()+1; ++j){
		// left boundary
		size_t ind_left = _grid.to_linear_point_index({0, j});
		ret[ind_left] = Vector(1, 0);
		// right boundary
		size_t ind_right = _grid.to_linear_point_index({_grid.nx(), j});
		double u_up, u_down;
		if (j < _grid.ny()){
			u_up = _u[_grid.yface_grid_index_i_jp(_grid.nx(), j)];
		} else {
			u_up = _u[_grid.yface_grid_index_i_jp(_grid.nx(), j-1)];
		}
		if (j > 0){
			u_down = _u[_grid.yface_grid_index_i_jp(_grid.nx(), j-1)];
		} else {
			u_down = _u[_grid.yface_grid_index_i_jp(_grid.nx(), j)];
		}
		double u = (u_up + u_down)/2.0;
		double v = _v[_grid.xface_grid_index_ip_j(_grid.nx()-1, j)];
		ret[ind_right] = Vector(u, v);
	}
	// internal
	for (size_t j=1; j<_grid.ny(); ++j)
	for (size_t i=1; i<_grid.nx(); ++i){
		size_t ind = _grid.to_linear_point_index({i, j});
		size_t ind_top = _grid.yface_grid_index_i_jp(i, j);
		size_t ind_bot = _grid.yface_grid_index_i_jp(i, j-1);
		size_t ind_left = _grid.xface_grid_index_ip_j(i-1, j);
		size_t ind_right = _grid.xface_grid_index_ip_j(i, j);
		double u, v;
		if (_grid.yface_type(ind_top) != RegularGrid2D::FaceType::Internal ||
		    _grid.yface_type(ind_bot) != RegularGrid2D::FaceType::Internal){
			u = 0.0;
		} else {
			u = 0.5*(_u[ind_top] + _u[ind_bot]);
		}
		if (_grid.xface_type(ind_left) != RegularGrid2D::FaceType::Internal ||
		    _grid.xface_type(ind_right) != RegularGrid2D::FaceType::Internal){
			v = 0.0;
		} else {
			v = 0.5*(_v[ind_left] + _v[ind_right]);
		}
		ret[ind] = Vector(u, v);
	}

	return ret;
}

ObstacleNonstat2DSimpleWorker::Coefficients ObstacleNonstat2DSimpleWorker::coefficients() const{
	double sum_cpx1 = 0;
	double sum_cpy1 = 0;
	double sum_cfx1 = 0;
	double sum_cfy1 = 0;
	double sum_nu1 = 0;
	double sum_cpx2 = 0;
	double sum_cpy2 = 0;
	double sum_cfx2 = 0;
	double sum_cfy2 = 0;
	double sum_nu2 = 0;

	// yfaces loop
	for (const RegularGrid2D::split_index_t& yface: _grid.boundary_yfaces()){
		if (yface[0] == 0){
			// input => ignore
		} else if (yface[0] == _grid.nx()){
			// output => ignore
		} else {
			double pnx, dvdn, dtdn;
			size_t left_cell = _grid.cell_centered_grid_index_ip_jp(yface[0]-1, yface[1]);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			int gamma_i = gamma_closest_to_cell(left_cell);
			double t_gamma = (gamma_i == 1) ? 0.5 : 1.0;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			size_t right_cell = _grid.cell_centered_grid_index_ip_jp(yface[0], yface[1]);

			if (_grid.is_active_cell(left_cell)){
				pnx = _p[left_cell];  
				dvdn = -v_ip_jp(yface[0]-1, yface[1]) / (_hx/2.0);
				dtdn = (t_gamma - _t[left_cell]) / (_hx / 2.0);
			} else if (_grid.is_active_cell(right_cell)){
				pnx = -_p[right_cell];  
				dvdn = -v_ip_jp(yface[0], yface[1]) / (_hx/2.0);
				dtdn = (t_gamma - _t[right_cell]) / (_hx / 2.0);
			} else {
				_THROW_UNREACHABLE_;
			}

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (gamma_i == 1) {
				sum_cpx1 += pnx * _hy;
				sum_cfy1 += dvdn * _hy;
				sum_nu1 += dtdn * _hy;
			}
			else {
				sum_cpx2 += pnx * _hy;
				sum_cfy2 += dvdn * _hy;
				sum_nu2 += dtdn * _hy;
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}

	// xfaces loop
	for (const RegularGrid2D::split_index_t& xface: _grid.boundary_xfaces()){
		if (xface[1] == 0){
			// bottom => ignore
		} else if (xface[1] == _grid.ny()){
			// top => ignore
		} else {
			double pny, dudn, dtdn;
			size_t bot_cell = _grid.cell_centered_grid_index_ip_jp(xface[0], xface[1]-1);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			int gamma_i = gamma_closest_to_cell(bot_cell);
			double t_gamma = (gamma_i == 1) ? 0.5 : 1.0;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			size_t top_cell = _grid.cell_centered_grid_index_ip_jp(xface[0], xface[1]);
			if (_grid.is_active_cell(bot_cell)){
				pny = _p[bot_cell];
				dudn = -u_ip_jp(xface[0], xface[1]-1)/(_hy/2.0);
				dtdn = (t_gamma - _t[bot_cell]) / (_hy / 2.0);
			} else if (_grid.is_active_cell(top_cell)){
				pny = -_p[top_cell];
				dudn = -u_ip_jp(xface[0], xface[1])/(_hy/2.0);
				dtdn = (t_gamma - _t[top_cell]) / (_hy / 2.0);
			} else {
				_THROW_UNREACHABLE_;
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (gamma_i == 1) {
				sum_cpx1 += pny * _hx;
				sum_cfy1 += dudn * _hx;
				sum_nu1 += dtdn * _hx;
			}
			else {
				sum_cpx2 += pny * _hx;
				sum_cfy2 += dudn * _hx;
				sum_nu2 += dtdn * _hx;
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}

	Coefficients coefs;
	coefs.cpx1 = 2.0*sum_cpx1;
	coefs.cpx2 = 2.0 * sum_cpx2;
	coefs.cpy1 = 2.0*sum_cpy1;
	coefs.cpy2 = 2.0 * sum_cpy2;
	coefs.cfx1 = -2.0/_Re*sum_cfx1;
	coefs.cfx2 = -2.0 / _Re * sum_cfx2;
	coefs.cfy1 = -2.0/_Re*sum_cfy1;
	coefs.cfy2 = -2.0 / _Re * sum_cfy2;
	coefs.cx1 = coefs.cpx1 + coefs.cfx1;
	coefs.cx2 = coefs.cpx2 + coefs.cfx2;
	coefs.cy1 = coefs.cpy1 + coefs.cfy1;
	coefs.cy2 = coefs.cpy2 + coefs.cfy2;
	coefs.nu1 = sum_nu1;
	coefs.nu2 = sum_nu2;
	return coefs;
}

namespace {

std::string convergence_report(double time, int it){
	std::ostringstream oss;
	oss << std::setprecision(1) << std::fixed;
	oss << "Time = " << std::setw(5) << time << " converged in " << it << " iterations" << std::endl;
	return oss.str();
}

}



//TEST_CASE("Obstacle 2D nonstationary, SIMPLE algorithm", "[obstacle2-nonstat-simple]"){
//	std::cout << std::endl << "--- cfd24_test [obstacle2-nonstat-simple] --- " << std::endl;
//
//	// problem parameters
//	double Re = 100;
//	double Pe = 100;
//	size_t n_unit = 10;  // partition per unit length
//	double time_step = 0.25;
//	double end_time = 5;
//	double E = 4.0;
//	size_t max_it = 10000;
//	double eps = 1e-0;
//
//	// worker initialization
//	RegularGrid2D grid(0, 12, -2, 2, 12*n_unit, 4*n_unit);
//	grid.deactivate_cells({2, -0.5}, {3, 0.5});
//	ObstacleNonstat2DSimpleWorker worker(Re, Pe, grid, E, time_step);
//	worker.initialize_saver(false, "obstacle2-nonstat");
//
//	// initial condition
//	worker.initialize();
//	worker.save_current_fields(0);
//
//	// iterations loop
//	for (double time=time_step; time<end_time+1e-6; time+=time_step){
//		size_t it = 0;
//		for (it=1; it < max_it; ++it){
//			double nrm = worker.step();
//
//			// break inner iterations if residual is low enough
//			if (nrm < eps){
//				break;
//			} else if (it == max_it -1) {
//				std::cout << "WARNING: internal SIMPLE interations not converged with nrm = "
//				          << nrm << std::endl;
//			}
//		}
//		// export solution each 1.0 time units
//		if (std::abs(time - round(time)) < 1e-6){
//			worker.save_current_fields(time);
//		}
//		std::cout << convergence_report(time, it);
//
//		// uvp_old = uvp
//		worker.to_next_time_step();
//	}
//
//	CHECK(worker.coefficients().cx == Approx(2.27).margin(1e-2));
//}

TEST_CASE("Obstacle 2D nonstationary, SIMPLE algorithm", "[obstacle2-nonstat-simple-2obstacles]") {
	std::cout << std::endl << "--- cfd24_test [obstacle2-nonstat-simple-2obstacles] --- " << std::endl;

	// problem parameters
	double Re = 100;
	double Pe = 100;
	size_t n_unit = 10;  // partition per unit length
	double time_step = 0.1;
	double end_time = 200;
	double E = 4.0;
	size_t max_it = 10000;
	double eps = 1e-1;

	// worker initialization
	RegularGrid2D grid(0, 12, -2, 2, 12 * n_unit, 4 * n_unit);
	grid.deactivate_cells({ 2, -0.7 }, { 2.5, 0.3 }); // first obstacle
	grid.deactivate_cells({ 4, -0.3 }, { 4.5, 0.7 }); // second obstacle
	ObstacleNonstat2DSimpleWorker worker(Re, Pe, grid, E, time_step);
	worker.initialize_saver(false, "obstacle2-nonstat-2");

	// initial condition
	worker.initialize();
	worker.save_current_fields(0);

	// iterations loop
	for (double time = time_step; time < end_time + 1e-6; time += time_step) {
		size_t it = 0;
		for (it = 1; it < max_it; ++it) {
			double nrm = worker.step();

			// break inner iterations if residual is low enough
			if (nrm < eps) {
				break;
			}
			else if (it == max_it - 1) {
				std::cout << "WARNING: internal SIMPLE interations not converged with nrm = "
					<< nrm << std::endl;
			}
		}
		// export solution each 1.0 time units
		if (std::abs(time - round(time)) < 1e-6) {
			worker.save_current_fields(time);
		}
		std::cout << convergence_report(time, it);

		// uvp_old = uvp
		worker.to_next_time_step();
	}

	//CHECK(worker.coefficients().cx == Approx(2.27).margin(1e-2));
}