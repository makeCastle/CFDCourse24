#include "cfd24_test.hpp"
#include "cfd24/grid/grid1d.hpp"
#include "cfd24/mat/csrmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/mat/lodmat.hpp"
#include "C:/git_repos/CFDCourse24/src/cfd24/grid/regular_grid2d.hpp"
//#include"C:/git_repos/CFDCourse24/src/cdf24/debug/printer.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Poisson 1D solver
///////////////////////////////////////////////////////////////////////////////

double kronecker_delta(size_t i, size_t j) {
	if (i == j)
	{
		return 1.0;
	}
	else {
		return 0.0;
	}
}


class TestPoisson1Worker{
public:
	// u(x) = sin(10*x^2)
	static double exact_solution(double x){
		return sin(10*x*x);
	}
	
	// -d^2 u(x)/d x^2
	static double exact_rhs(double x){
		return 400*x*x*sin(10*x*x) - 20*cos(10*x*x);
	}

	TestPoisson1Worker(size_t n_cells): grid(0, 1, n_cells){}

	// returns norm2(u - u_exact)
	double solve(){
		// 1. build SLAE
		CsrMatrix mat = approximate_lhs();
		std::vector<double> rhs = approximate_rhs();

		// 2. solve SLAE
		AmgcMatrixSolver solver;
		solver.set_matrix(mat);
		solver.solve(rhs, u);

		// 3. compute norm2
		return compute_norm2();
	}

	// saves numerical and exact solution into the vtk format
	void save_vtk(const std::string& filename){
		// save grid
		grid.save_vtk(filename);
		
		// save numerical solution
		VtkUtils::add_point_data(u, "numerical", filename);

		// save exact solution
		std::vector<double> exact(grid.n_points());
		for (size_t i=0; i<grid.n_points(); ++i){
			exact[i] = exact_solution(grid.point(i).x());
		}
		VtkUtils::add_point_data(exact, "exact", filename);
	}

private:
	const Grid1D grid;
	std::vector<double> u;

	CsrMatrix approximate_lhs() const{
		// constant h = x[1] - x[0]
		double h = grid.point(1).x() - grid.point(0).x();
		
		// fill using 'easy-to-construct' sparse matrix format
		LodMatrix mat(grid.n_points());
		mat.add_value(0, 0, 1);
		mat.add_value(grid.n_points()-1, grid.n_points()-1, 1);
		double diag = 2.0/h/h;
		double nondiag = -1.0/h/h;
		for (size_t i=1; i<grid.n_points()-1; ++i){
			mat.add_value(i, i-1, nondiag);
			mat.add_value(i, i+1, nondiag);
			mat.add_value(i, i, diag);
		}

		// return 'easy-to-use' sparse matrix format
		return mat.to_csr();
	}

	std::vector<double> approximate_rhs() const{
		std::vector<double> ret(grid.n_points());
		ret[0] = exact_solution(grid.point(0).x());
		ret[grid.n_points()-1] = exact_solution(grid.point(grid.n_points()-1).x());
		for (size_t i=1; i<grid.n_points()-1; ++i){
			ret[i] = exact_rhs(grid.point(i).x());
		}
		return ret;
	}

	double compute_norm2() const{
		// weights
		double h = grid.point(1).x() - grid.point(0).x();
		std::vector<double> w(grid.n_points(), h);
		w[0] = w[grid.n_points()-1] = h/2;

		// sum
		double sum = 0;
		for (size_t i=0; i<grid.n_points(); ++i){
			double diff = u[i] - exact_solution(grid.point(i).x());
			sum += w[i]*diff*diff;
		}

		double len = grid.point(grid.n_points()-1).x() - grid.point(0).x();
		return std::sqrt(sum / len);
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class TestPoisson2Worker {
public:
	// u = cos(10 * x * y) + sin(10 * x * y)
	static double exact_solution(double x, double y) {
		return cos(10 * x * y) + sin(10 * x * y);
	}

	static double exact_rhs(double x, double y) {
		return  100 * (x * x * cos(10 * y * x) + y * y * cos(10 * y * x) + x * x * sin(10 * y * x) + y * y * sin(10 * y * x));
	}

	TestPoisson2Worker(size_t n_cellsX, size_t n_cellsY) : _grid(0, 1, 0, 1, n_cellsX, n_cellsY), n_cellsX(n_cellsX), n_cellsY(n_cellsY), _hx(1.0/n_cellsX), _hy(1.0/n_cellsY) {}

	double solve() {
		// 1. build SLAE
		CsrMatrix mat = approximate_lhs();
		std::vector<double> rhs = approximate_rhs();

		// 2. solve SLAE
		AmgcMatrixSolver solver;
		solver.set_matrix(mat);
		solver.solve(rhs, _u);

		// 3. compute norm2
		return compute_norm2();
	}

	// saves numerical and exact solution into the vtk format
	void save_vtk(const std::string& filename) {
		// save grid
		_grid.save_vtk(filename);

		// save numerical solution
		VtkUtils::add_point_data(_u, "numerical", filename);

		// save exact solution
		std::vector<double> exact(_grid.n_points());
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			exact[i] = exact_solution(_grid.point(i).x(), _grid.point(i).y());
		}
		VtkUtils::add_point_data(exact, "exact", filename);
	}



private:
	const RegularGrid2D _grid;
	std::vector<double> _u;
	size_t n_cellsX;
	size_t n_cellsY;
	double _hx;
	double _hy;

	CsrMatrix approximate_lhs() const {
		// fill using 'easy-to-construct' sparse matrix format
		LodMatrix mat(_grid.n_points());

		for (size_t i_point = 0; i_point < _grid.n_points(); i_point++) {
			// i_point
			RegularGrid2D::split_index_t ij = _grid.to_split_point_index(i_point);
			size_t i = ij[0];
			size_t j = ij[1];
			//��� �������
			if (i == 0 || i == n_cellsX || j == n_cellsY || j == 0) {
				mat.add_value(i_point, i_point, 1);
			}


			else {
				//k[i-1,j]
				RegularGrid2D::split_index_t ij_im_j{ i - 1, j };
				size_t k_im_j = _grid.to_linear_point_index(ij_im_j);
				mat.add_value(i_point, k_im_j, -1.0 / _hx / _hx);
				//k[i,j]
				RegularGrid2D::split_index_t ij_i_j{ i, j };
				size_t k_i_j = _grid.to_linear_point_index(ij_i_j);
				mat.add_value(i_point, k_i_j, 2.0 / _hx / _hx + 2.0 / _hy / _hy);
				//k[i+1,j]
				RegularGrid2D::split_index_t ij_ip_j{ i + 1, j };
				size_t k_ip_j = _grid.to_linear_point_index(ij_ip_j);
				mat.add_value(i_point, k_ip_j, -1.0 / _hx / _hx);
				//k[i,j-1]
				RegularGrid2D::split_index_t ij_i_jm{ i , j - 1 };
				size_t k_i_jm = _grid.to_linear_point_index(ij_i_jm);
				mat.add_value(i_point, k_i_jm, -1.0 / _hy / _hy);
				//k[i,j+1]
				RegularGrid2D::split_index_t ij_i_jp{ i, j + 1 };
				size_t k_i_jp = _grid.to_linear_point_index(ij_i_jp);
				mat.add_value(i_point, k_i_jp, -1.0 / _hy / _hy);
			}

		}
		// return 'easy-to-use' sparse matrix format
		return mat.to_csr();
	}

	std::vector<double> approximate_rhs() const {
		std::vector<double> ret(_grid.n_points());
		for (size_t i_point = 0; i_point < _grid.n_points(); i_point++) {
			Point p = _grid.point(i_point);
			// i_point
			RegularGrid2D::split_index_t ij = _grid.to_split_point_index(i_point);
			size_t i = ij[0];
			size_t j = ij[1];
			//
			if (i == 0 || i == n_cellsX || j == n_cellsY || j == 0) {
				ret[i_point] = exact_solution(p.x(), p.y());
			}


			else {
				ret[i_point] = exact_rhs(p.x(), p.y());
			}

		}
		return ret;
	}

	double compute_norm2() const {
		// sum
		double sum = 0;
		for (size_t i = 0; i < n_cellsX + 1; ++i) {
			for (size_t j = 0; j < n_cellsY + 1; ++j)
			{
				double Ek = _hx * _hy;
				size_t k = _grid.to_linear_point_index({ i,j });
				if ((i == 0 && j == 0) || (i == n_cellsX && j == 0) || (i == 0 && j == n_cellsY) || (i == n_cellsX && j == n_cellsY)) {
					Ek = Ek / 4;
				}
				else if (i == 0 || j == 0 || i == n_cellsX || j == n_cellsY) {
					Ek = Ek / 2;
				}
				double diff = _u[k] - exact_solution(_grid.point(k).x(), _grid.point(k).y());
				sum += Ek * diff * diff;
			}
		}
		return sqrt(sum);
	}
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Poisson 1D solver", "[poisson1]"){
	std::cout << std::endl << "--- cfd24_test [poisson1] --- " << std::endl;

	// precalculated norm2 results for some n_cells values
	// used for CHECK procedures
	std::map<size_t, double> norm2_for_compare{
		{10, 0.179124},
		{100, 0.00158055},
		{1000, 1.57849e-05},
	};

	// loop over n_cells value
	for (size_t n_cells: {10, 20, 50, 100, 200, 500, 1000}){
		// build test solver
		TestPoisson1Worker worker(n_cells);

		// solve and find norm2
		double n2 = worker.solve();

		// save into poisson1_ncells={n_cells}.vtk
		worker.save_vtk("poisson1_ncells=" + std::to_string(n_cells) + ".vtk");

		// print (N_CELLS, NORM2) table entry
		std::cout << n_cells << " " << n2 << std::endl;

		// CHECK if result for this n_cells
		// presents in the norm2_for_compare dictionary
		auto found = norm2_for_compare.find(n_cells);
		if (found != norm2_for_compare.end()){
			CHECK(n2 == Approx(found->second).margin(1e-6));
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Poisson 2D solver", "[poisson2]") {
	std::cout << std::endl << "--- cfd24_test [poisson2] --- " << std::endl;

	for (size_t n_cells : {10, 20, 50, 100, 200, 500}) {
		// build test solver
		TestPoisson2Worker worker(n_cells, n_cells);


		// solve and find norm2
		double n2 = worker.solve();

		// save into poisson1_ncells={n_cells}.vtk
		worker.save_vtk("poisson2_ncells=" + std::to_string(n_cells) + ".vtk");

		// print (N_CELLS, NORM2) table entry
		std::cout << n_cells << " " << n2 << std::endl;

		// CHECK if result for this n_cells
		// presents in the norm2_for_compare dictionary
		/*auto found = norm2_for_compare.find(n_cells);
		if (found != norm2_for_compare.end()){
			CHECK(n2 == Approx(found->second).margin(1e-6));
		}*/

	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////