#include "cfd24_test.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/mat/lodmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"

using namespace cfd;


struct Cavern2DSimpleWorker{
	Cavern2DSimpleWorker(double Re, size_t n_cells, double alpha_p):
		_grid(0, 1, 0, 1, n_cells, n_cells),
		_cc_grid(_grid.cell_centered_grid()),
		_xe_grid(_grid.xedge_centered_grid()),
		_ye_grid(_grid.yedge_centered_grid()),
		_hx(1.0/n_cells),
		_hy(1.0/n_cells),
		_Re(Re),
		_alpha_p(alpha_p),
		_p(_cc_grid.n_points(), 0.0),
		_u(_ye_grid.n_points(), 0.0),
		_v(_xe_grid.n_points(), 0.0),
		_du(_ye_grid.n_points(), 0.0),
		_dv(_xe_grid.n_points(), 0.0)
	{
	}

	double step(double tau){
		_current_time += tau;

		// U-star
		std::vector<double> u_star = compute_u_star(tau);
		std::vector<double> v_star = compute_v_star(tau);
		// Pressure correction
		std::vector<double> p_stroke = compute_p_stroke(tau, u_star, v_star);
		// Velocity correction
		std::vector<double> u_stroke = compute_u_stroke(tau, p_stroke);
		std::vector<double> v_stroke = compute_v_stroke(tau, p_stroke);
		// Set final values
		double nrm_u = set_final_u(u_star, u_stroke);
		double nrm_v = set_final_v(v_star, v_stroke);
		set_final_p(p_stroke);

		return std::max(nrm_u, nrm_v) / tau;
	}

	void save_vtk(std::string file_u, std::string file_v, std::string file_p){
		// pressure
		_cc_grid.save_vtk(file_p);
		VtkUtils::add_point_data(_p, "pressure", file_p);
		// u
		_ye_grid.save_vtk(file_u);
		VtkUtils::add_point_data(_u, "velocity-x", file_u);
		// v
		_xe_grid.save_vtk(file_v);
		VtkUtils::add_point_data(_v, "velocity-y", file_v);
	}

	double current_time() const{
		return _current_time;
	}

private:
	const RegularGrid2D _grid;
	const RegularGrid2D _cc_grid;
	const RegularGrid2D _xe_grid;
	const RegularGrid2D _ye_grid;
	const double _hx;
	const double _hy;
	const double _Re;
	const double _alpha_p;
	double _current_time = 0;

	std::vector<double> _p;
	std::vector<double> _u;
	std::vector<double> _v;
	std::vector<double> _du;
	std::vector<double> _dv;
	
	std::vector<double> compute_u_star(double tau){
		std::vector<double> rhs(_u.size(), 0);
		LodMatrix mat(_u.size());
		for (size_t ipoint=0; ipoint<_grid.n_points(); ++ipoint){
			std::array<size_t, 2> ij = _grid.to_split_point_index(ipoint);
			size_t i  = ij[0];
			size_t j  = ij[1];
			if (j == _grid.ny()){
				// out of area
				continue;
			}
			size_t row_index = _grid.yedge_grid_index_i_jp(i, j);   //[i, j+1/2] linear index in u grid
			bool is_left = (i == 0);
			bool is_right = (i == _grid.nx());
			bool is_bottom = (j == 0);
			bool is_top = (j == _grid.ny() - 1);
			if (is_left || is_right){
				// === left/right boundary condition: dirichlet value
				mat.set_value(row_index, row_index, 1.0);
				rhs[row_index] = 0.0;
			} else {
				// === internal
				size_t x_plus_index  = _grid.yedge_grid_index_i_jp(i+1, j); //(i+1, j+1/2) index
				size_t x_minus_index = _grid.yedge_grid_index_i_jp(i-1, j); //(i-1, j+1/2) index
				size_t y_plus_index = (size_t)(-1);
				if (!is_top){
					y_plus_index  = _grid.yedge_grid_index_i_jp(i, j+1); //(i, j+3/2) index
				}
				size_t y_minus_index = (size_t)(-1);
				if (!is_bottom){
					y_minus_index = _grid.yedge_grid_index_i_jp(i, j-1); //(i, j-1/2) index
				}
				// u_i_jp
				mat.set_value(row_index, row_index, 1.0);
				// + tau * d (u_prev*u) / dx
				double u_prev_plus   = u_ip_jp(i, j);  //_u[i+1/2, j+1/2]
				double u_prev_minus  = u_ip_jp(i-1, j);//_u[i-1/2, j+1/2]
				mat.add_value(row_index, row_index,
					tau/2.0/_hx*(u_prev_plus - u_prev_minus));
				mat.add_value(row_index, x_plus_index,
					tau/2.0/_hx*(u_prev_plus));
				mat.add_value(row_index, x_minus_index,
					tau/2.0/_hx*(-u_prev_minus));
				// + tau * d (v_prev*u) / dy
				double v_prev_plus   = v_i_j(i, j+1); // _v[i,j+1]
				double v_prev_minus  = v_i_j(i, j);   // _v[i,j]
				mat.add_value(row_index, row_index,
					tau/2.0/_hy*(v_prev_plus - v_prev_minus));
				if (is_top){
					mat.add_value(row_index, row_index,
						-tau/2.0/_hy*(v_prev_plus));
					rhs[row_index] -= 2.0*tau/2.0/_hy*(v_prev_plus);
				} else {
					mat.add_value(row_index, y_plus_index,
						tau/2.0/_hy*(v_prev_plus));
				}
				if (is_bottom){
					mat.add_value(row_index, row_index,
						-tau/2.0/_hy*(-v_prev_minus));
				} else {
					mat.add_value(row_index, y_minus_index,
						tau/2.0/_hy*(-v_prev_minus));
				}
				// - tau/Re[d^2 u/d x^2]
				mat.add_value(row_index, row_index, 2.0*tau/_Re/_hx/_hx);
				mat.add_value(row_index, x_plus_index, -tau/_Re/_hx/_hx);
				mat.add_value(row_index, x_minus_index, -tau/_Re/_hx/_hx);
				// - tau/Re[d^2 u/d y^2]
				mat.add_value(row_index, row_index, 2.0*tau/_Re/_hy/_hy);
				if (is_top){
					mat.add_value(row_index, row_index, tau/_Re/_hy/_hy);
					rhs[row_index] += 2*tau/_Re/_hy/_hy;
				} else {
					mat.add_value(row_index, y_plus_index, -tau/_Re/_hy/_hy);
				}
				if (is_bottom){
					mat.add_value(row_index, row_index, tau/_Re/_hy/_hy);
				} else {
					mat.add_value(row_index, y_minus_index, -tau/_Re/_hy/_hy);
				}
				// = u_prev
				rhs[row_index] += _u[row_index];
				//   -tau * dp/dx
				rhs[row_index] -= tau/_hx*(p_ip_jp(i, j) - p_ip_jp(i-1, j));
			}
		}

		CsrMatrix csr = mat.to_csr();
		AmgcMatrixSolver slv;
		slv.set_matrix(csr);
		std::vector<double> u_star(_u);
		slv.solve(rhs, u_star);

		// fill du
		_du = csr.diagonal();
		for (double& v: _du) v = 1.0/v;

		return u_star;
	}

	std::vector<double> compute_v_star(double tau){
		std::vector<double> rhs(_v.size(), 0);
		LodMatrix mat(_v.size());
		for (size_t ipoint=0; ipoint<_grid.n_points(); ++ipoint){
			std::array<size_t, 2> ij = _grid.to_split_point_index(ipoint);
			size_t i  = ij[0];
			size_t j  = ij[1];
			if (i == _grid.nx()){
				// out of area
				continue;
			}
			size_t row_index = _grid.xedge_grid_index_ip_j(i, j);   //[i+1/2, j] linear index in v grid
			bool is_left = (i == 0);
			bool is_right = (i == _grid.nx()-1);
			bool is_bottom = (j == 0);
			bool is_top = (j == _grid.ny());
			if (is_bottom || is_top){
				// === top/bottom boundary condition: dirichlet value = 0
				mat.set_value(row_index, row_index, 1.0);
				rhs[row_index] = 0.0;
			} else {
				// === internal
				size_t x_plus_index  = (size_t)(-1);
				if (!is_right){
					x_plus_index = _grid.xedge_grid_index_ip_j(i+1, j); //(i+3/2, j) index
				}
				size_t x_minus_index = (size_t)(-1);
				if (!is_left){
					x_minus_index =_grid.xedge_grid_index_ip_j(i-1, j); //(i-1/2, j) index
				}
				size_t y_plus_index = _grid.xedge_grid_index_ip_j(i, j+1);  //(i+1/2, j+1) index
				size_t y_minus_index = _grid.xedge_grid_index_ip_j(i, j-1);  //(i+1/2, j-1) index
				// v_ip_j
				mat.set_value(row_index, row_index, 1.0);
				// + tau * d (u_prev*v) / dx
				double u_prev_plus   = u_i_j(i+1, j);  //_u[i+1, j]
				double u_prev_minus  = u_i_j(i, j);    //_u[i, j]
				mat.add_value(row_index, row_index,
					tau/2.0/_hx*(u_prev_plus - u_prev_minus));
				if (is_right){
					mat.add_value(row_index, row_index,
						-tau/2.0/_hx*(u_prev_plus));
				} else {
					mat.add_value(row_index, x_plus_index,
						tau/2.0/_hx*(u_prev_plus));
				}
				if (is_left){
					mat.add_value(row_index, row_index,
						-tau/2.0/_hx*(-u_prev_minus));
				} else {
					mat.add_value(row_index, x_minus_index,
						tau/2.0/_hx*(-u_prev_minus));
				}
				// + tau * d (v_prev*v) / dy
				double v_prev_plus   = v_ip_jp(i, j);   // _v[i+1/2, j+1/2]
				double v_prev_minus  = v_ip_jp(i, j-1); // _v[i+1/2, j-1/2]
				mat.add_value(row_index, row_index,
					tau/2.0/_hy*(v_prev_plus - v_prev_minus));
				mat.add_value(row_index, y_plus_index,
					tau/2.0/_hy*(v_prev_plus));
				mat.add_value(row_index, y_minus_index,
					tau/2.0/_hy*(-v_prev_minus));
				// - tau/Re[d^2 v/d x^2]
				mat.add_value(row_index, row_index, 2.0*tau/_Re/_hx/_hx);
				if (is_right){
					mat.add_value(row_index, row_index, tau/_Re/_hx/_hx);
				} else {
					mat.add_value(row_index, x_plus_index, -tau/_Re/_hx/_hx);
				}
				if (is_left){
					mat.add_value(row_index, row_index, tau/_Re/_hx/_hx);
				} else {
					mat.add_value(row_index, x_minus_index, -tau/_Re/_hx/_hx);
				}
				// - tau/Re[d^2 v/d y^2]
				mat.add_value(row_index, row_index, 2.0*tau/_Re/_hy/_hy);
				mat.add_value(row_index, y_plus_index, -tau/_Re/_hy/_hy);
				mat.add_value(row_index, y_minus_index, -tau/_Re/_hy/_hy);
				// = v_prev
				rhs[row_index] += _v[row_index];
				//   -tau * dp/dy
				rhs[row_index] -= tau/_hy*(p_ip_jp(i, j) - p_ip_jp(i, j-1));
			}
		}

		CsrMatrix csr = mat.to_csr();
		AmgcMatrixSolver slv;
		slv.set_matrix(csr);
		std::vector<double> v_star(_v);
		slv.solve(rhs, v_star);

		// fill dv
		_dv = csr.diagonal();
		for (double& v: _dv) v = 1.0/v;

		return v_star;
	}

	std::vector<double> compute_p_stroke(double tau, const std::vector<double>& u_star, const std::vector<double>& v_star){
		for (double& v: _du) v = 1.0;
		for (double& v: _dv) v = 1.0;

		LodMatrix mat(_p.size());
		std::vector<double> rhs(_p.size(), 0.0);
		for (size_t ipoint = 0; ipoint < _grid.n_points(); ++ipoint){
			std::array<size_t, 2> ij = _grid.to_split_point_index(ipoint);
			size_t i  = ij[0];
			size_t j  = ij[1];
			if (i == _grid.nx() || j == _grid.ny()){
				continue;
			}
			bool is_left = (i==0);
			bool is_right = (i==_grid.nx()-1);
			bool is_bottom = (j==0);
			bool is_top = (j==_grid.ny() - 1);

			size_t ind0 = _grid.cell_centered_grid_index_ip_jp(i, j);
			// x
			if (!is_right){
				double coef = tau*_du[_grid.yedge_grid_index_i_jp(i+1, j)]/_hx/_hx;
				size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i+1, j);
				mat.add_value(ind0, ind0, coef);
				mat.add_value(ind0, ind1, -coef);
			}
			if (!is_left){
				double coef = tau*_du[_grid.yedge_grid_index_i_jp(i, j)]/_hx/_hx;
				size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i-1, j);
				mat.add_value(ind0, ind0, coef);
				mat.add_value(ind0, ind1, -coef);
			}
			// y
			if (!is_top){
				double coef = tau*_dv[_grid.xedge_grid_index_ip_j(i, j+1)]/_hy/_hy;
				size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i, j+1);
				mat.add_value(ind0, ind0, coef);
				mat.add_value(ind0, ind1, -coef);
			}
			if (!is_bottom){
				double coef = tau*_dv[_grid.xedge_grid_index_ip_j(i, j)]/_hy/_hy;
				size_t ind1 = _grid.cell_centered_grid_index_ip_jp(i, j-1);
				mat.add_value(ind0, ind0, coef);
				mat.add_value(ind0, ind1, -coef);
			}
			// rhs
			size_t ind_left = _grid.yedge_grid_index_i_jp(i, j);
			size_t ind_right = _grid.yedge_grid_index_i_jp(i+1, j);
			size_t ind_bot = _grid.xedge_grid_index_ip_j(i, j);
			size_t ind_top = _grid.xedge_grid_index_ip_j(i, j+1);
			rhs[ind0] = -(u_star[ind_right] - u_star[ind_left])/_hx - (v_star[ind_top] - v_star[ind_bot])/_hy;
		}

		CsrMatrix csr = mat.to_csr();
		AmgcMatrixSolver slv;
		slv.set_matrix(csr);
		std::vector<double> p_stroke;
		slv.solve(rhs, p_stroke);
		return p_stroke;
	}

	std::vector<double> compute_u_stroke(double tau, const std::vector<double>& p_stroke){
		std::vector<double> u_stroke(_u.size(), 0.0);
		for (size_t ipoint = 0; ipoint < _grid.n_points(); ++ipoint){
			std::array<size_t, 2> ij = _grid.to_split_point_index(ipoint);
			size_t i  = ij[0];
			size_t j  = ij[1];
			if (j == _grid.ny()){
				continue;
			}
			bool is_left = (i == 0);
			bool is_right = (i == _grid.nx());
			size_t ind0 = _grid.yedge_grid_index_i_jp(i, j);
			if (is_left || is_right){
				u_stroke[ind0] = 0;
			} else {
				size_t ind_plus  = _grid.cell_centered_grid_index_ip_jp(i, j);
				size_t ind_minus = _grid.cell_centered_grid_index_ip_jp(i-1, j);
				u_stroke[ind0] = -tau * _du[ind0] * (p_stroke[ind_plus] - p_stroke[ind_minus])/_hx;
			}
		}
		return u_stroke;
	}

	std::vector<double> compute_v_stroke(double tau, const std::vector<double>& p_stroke){
		std::vector<double> v_stroke(_v.size(), 0.0);
		for (size_t ipoint = 0; ipoint < _grid.n_points(); ++ipoint){
			std::array<size_t, 2> ij = _grid.to_split_point_index(ipoint);
			size_t i  = ij[0];
			size_t j  = ij[1];
			if (i == _grid.nx()){
				continue;
			}
			bool is_bottom = (j == 0);
			bool is_top = (j == _grid.ny());
			size_t ind0 = _grid.xedge_grid_index_ip_j(i, j);
			if (is_bottom || is_top){
				v_stroke[ind0] = 0;
			} else {
				size_t ind_plus  = _grid.cell_centered_grid_index_ip_jp(i, j);
				size_t ind_minus = _grid.cell_centered_grid_index_ip_jp(i, j-1);
				v_stroke[ind0] = -tau * _dv[ind0] * (p_stroke[ind_plus] - p_stroke[ind_minus])/_hy;
			}
		}
		return v_stroke;
	}

	double set_final_u(const std::vector<double>& u_star, const std::vector<double>& u_stroke){
		double nrm = 0;
		for (size_t i=0; i<_u.size(); ++i){
			nrm = std::max(nrm, std::abs(_u[i] - u_star[i] - u_stroke[i]));
			_u[i] = u_star[i] + u_stroke[i];
		}
		return nrm;
	}

	double set_final_v(const std::vector<double>& v_star, const std::vector<double>& v_stroke){
		double nrm = 0;
		for (size_t i=0; i<_v.size(); ++i){
			nrm = std::max(nrm, std::abs(_v[i] - v_star[i] - v_stroke[i]));
			_v[i] = v_star[i] + v_stroke[i];
		}
		return nrm;
	}

	double set_final_p(const std::vector<double>& p_stroke){
		double nrm = 0;
		for (size_t i=0; i<_p.size(); ++i){
			nrm = std::max(nrm, std::abs(p_stroke[i]));
			_p[i] = _p[i] + _alpha_p * p_stroke[i];
		}
		return nrm;
	}

	double u_i_j(size_t i, size_t j) const{
		double ind0 = _grid.yedge_grid_index_i_jp(i, j);
		double ind1 = _grid.yedge_grid_index_i_jp(i, j-1);
		return (_u[ind0] + _u[ind1])/2.0;
	}

	double u_ip_jp(size_t i, size_t j) const{
		double ind0 = _grid.yedge_grid_index_i_jp(i, j);
		double ind1 = _grid.yedge_grid_index_i_jp(i+1, j);
		return (_u[ind0] + _u[ind1])/2.0;
	}

	double v_i_j(size_t i, size_t j) const{
		double ind0 = _grid.xedge_grid_index_ip_j(i, j);
		double ind1 = _grid.xedge_grid_index_ip_j(i-1, j);
		return (_v[ind0] + _v[ind1])/2.0;
	}
	double v_ip_jp(size_t i, size_t j) const{
		double ind0 = _grid.xedge_grid_index_ip_j(i, j);
		double ind1 = _grid.xedge_grid_index_ip_j(i, j+1);
		return (_v[ind0] + _v[ind1])/2.0;
	}

	double p_ip_jp(size_t i, size_t j) const{
		return _p[_grid.cell_centered_grid_index_ip_jp(i, j)];
	}
};

TEST_CASE("Cavern 2D, SIMPLE algorithm", "[cavern2-simple]"){
	std::cout << std::endl << "--- cfd24_test [cavern2-simple] --- " << std::endl;

	double Re = 100;
	double tau = 1.0;
	double alpha = 0.8;
	size_t n_cells = 30;
	size_t max_it = 10000;
	double eps = 1e-2;

	// worker initialization
	Cavern2DSimpleWorker worker(Re, n_cells, alpha);

	// saver initialization
	VtkUtils::TimeDependentWriter writer_u("cavern2d-u");
	VtkUtils::TimeDependentWriter writer_v("cavern2d-v");
	VtkUtils::TimeDependentWriter writer_p("cavern2d-p");
	std::string out_filename_u = writer_u.add(0);
	std::string out_filename_v = writer_v.add(0);
	std::string out_filename_p = writer_p.add(0);
	worker.save_vtk(out_filename_u, out_filename_v, out_filename_p);

	size_t it = 0;
	for (it=0; it < max_it; ++it){
		double nrm = worker.step(tau);

		// print norm
		std::cout << it << " " << nrm << std::endl;

		// export solution to vtk
		out_filename_u = writer_u.add(worker.current_time());
		out_filename_v = writer_v.add(worker.current_time());
		out_filename_p = writer_p.add(worker.current_time());
		worker.save_vtk(out_filename_u, out_filename_v, out_filename_p);

		if (nrm < eps){
			break;
		}
	}
	CHECK(it == 10);
}
