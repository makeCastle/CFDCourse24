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
	// u(x) = sin(10*x^2)
	static double exact_solution(double x, double y) {
		return sin(10 * x * y);
	}

	// -d^2 u(x)/d x^2
	static double exact_rhsX(double x, double y) {
		return -100 * y * y * sin(10 * x * y);
	}

	// -d^2 u(x)/d y^2
	static double exact_rhsY(double x, double y) {
		return -100 * x * x * sin(10 * x * y);
	}

	TestPoisson2Worker(size_t n_cellsX, size_t n_cellsY) : grid(0.0, 1.0, 0.0, 1.0, n_cellsX, n_cellsY), _hx(1.0 / n_cellsX), _hy(1.0 / n_cellsY), _nx(n_cellsX + 1), _ny(n_cellsY + 1) {}

	// returns norm2(u - u_exact)
	double solve() {
		// 1. build SLAE
		CsrMatrix mat = approximate_lhs2();
		//dbg::print(mat);

		std::vector<double> rhs = approximate_rhs2();

		// 2. solve SLAE
		AmgcMatrixSolver solver;
		solver.set_matrix(mat);
		solver.solve(rhs, u);

		// 3. compute norm2
		return compute_norm2D();
	}

	// saves numerical and exact solution into the vtk format
	void save_vtk(const std::string& filename) {
		// save grid
		grid.save_vtk(filename);

		// save numerical solution
		VtkUtils::add_point_data(u, "numerical", filename);

		// save exact solution
		std::vector<double> exact(grid.n_points());
		for (size_t i = 0; i < grid.n_points(); ++i) {
			exact[i] = exact_solution(grid.point(i).x(), grid.point(i).y());
		}
		VtkUtils::add_point_data(exact, "exact", filename);
	}

protected:
	double _hx;
	double _hy;
	size_t _ny;
	size_t _nx;
private:
	const RegularGrid2D grid;
	std::vector<double> u;

	CsrMatrix approximate_lhs2() const {
		// constant h = x[1] - x[0]
		//size_t a = grid.point(grid.n_points() - 1).x() / (grid.point(1).x() - grid.point(0).x()); // number of last point in first row
		//size_t b = grid.point(grid.n_points() - 1).y() / (grid.point(2 * a + 1).y() - grid.point(a).y()); // number of last point in first collumn

		//double hx = grid.point(1).x() - grid.point(0).x();
		//double hy = grid.point(2 * a + 1).y() - grid.point(a).y();

		// fill using 'easy-to-construct' sparse matrix format
		LodMatrix mat(grid.n_points());
		//mat.add_value(0, 0, 1);
		//mat.add_value(grid.n_points() - 1, grid.n_points() - 1, 1);

		double nondiag1 = -1.0 / (_hx * _hx);
		double nondiag2 = -1.0 / (_hy * _hy);
		double diag = 2.0 / (_hx * _hx) + 2.0 / (_hy * _hy);
		double nondiag3 = -1.0 / (_hy * _hy);
		double nondiag4 = -1.0 / (_hx * _hx);
		
		/*mat.add_value(1, 0, diag); mat.add_value(1, 1, nondiag3); mat.add_value(1, 2, nondiag4);
		mat.add_value(grid.n_points() - 2, grid.n_points() - 3, nondiag1); mat.add_value(grid.n_points() - 2, grid.n_points() - 2, nondiag2); mat.add_value(grid.n_points() - 2, grid.n_points() - 1, diag);
		mat.add_value(2, 0, nondiag2); mat.add_value(2, 1, diag); mat.add_value(2, 2, nondiag3); mat.add_value(2, 3, nondiag4);
		mat.add_value(grid.n_points() - 3, grid.n_points() - 4, nondiag1); mat.add_value(grid.n_points() - 3, grid.n_points() - 3, nondiag2); 
		mat.add_value(grid.n_points() - 3, grid.n_points() - 2, diag); mat.add_value(grid.n_points() - 3, grid.n_points() - 1, nondiag3);*/

		//mat.add_value(1, 0, nondiag2); mat.add_value(1, 1, diag); mat.add_value(1, 2, nondiag3); mat.add_value(1, 3, nondiag4);
		//mat.add_value(grid.n_points() - 2, grid.n_points() - 4, nondiag1); mat.add_value(grid.n_points() - 2, grid.n_points() - 3, nondiag2);
		//mat.add_value(grid.n_points() - 2, grid.n_points() - 2, diag); mat.add_value(grid.n_points() - 2, grid.n_points() - 1, nondiag3);
		
		mat.add_value(0, 0, 1);
		mat.add_value(grid.n_points() - 1, grid.n_points() - 1, 1);
		for (size_t i = 1; i < _ny - 1; ++i)
		{
			for (size_t j = 1; j < _nx - 1; ++j)
			{
				size_t index = grid.to_linear_point_index({ j, i });
				mat.add_value(index, index - 1, nondiag2);
				mat.add_value(index, index + 1, nondiag3);
				mat.add_value(index, index, diag);
				mat.add_value(index - 1, index, nondiag1);
				mat.add_value(index + 1, index, nondiag4);
			}
		}

		//for (size_t i = 0; i < grid.n_points(); ++i) {
		//	cfd::RegularGrid2D::split_index_t irow_col = grid.to_split_point_index(i);
		//	size_t irow = irow_col[1]; size_t icol = irow_col[0];
		//	if (i == 24)
		//	{
		//		int a = 0;
		//	}
		//	//// BEGIN boundary conditions
		//	//if (irow == 0)
		//	//{
		//	//	for (size_t j = 0; j < _nx; ++j)
		//	//	{
		//	//		mat.add_value(i, j * _nx, kronecker_delta(irow, icol));
		//	//	}
		//	//}
		//	//else if (icol == 0 || icol == _nx - 1)
		//	//{
		//	//	for (size_t j = 0; j < _ny; ++j)
		//	//	{
		//	//		mat.add_value(i + j, i, kronecker_delta(irow, icol));
		//	//	}
		//	//}
		//	//// END boundary conditions
		//	
		//	if ((i == 0) || (i == grid.n_points() - 1))
		//	{
		//		mat.add_value(i, i, 1);
		//	}

		//	if ((irow > 0) && (irow < _ny - 1))
		//	{
		//		if (irow == 1)
		//		{
		//			mat.add_value(i, i - 1, nondiag2);
		//			mat.add_value(i, i, diag);
		//			mat.add_value(i, i + 1, nondiag3);
		//			mat.add_value(i, i + _nx - 1, nondiag4);
		//		}
		//		else if (irow == _ny - 2)
		//		{
		//			mat.add_value(i, i - _nx + 1, nondiag1);
		//			mat.add_value(i, i - 1, nondiag2);
		//			mat.add_value(i, i, diag);
		//			mat.add_value(i, i + 1, nondiag3);
		//		}
		//		else {
		//			mat.add_value(i, i - _nx + 1, nondiag1);
		//			mat.add_value(i, i - 1, nondiag2);
		//			mat.add_value(i, i, diag);
		//			mat.add_value(i, i + 1, nondiag3);
		//			mat.add_value(i, i + _nx - 1, nondiag4);
		//		}
		//	}
		//	

		//	

		//}

		dbg::print(mat);
		// return 'easy-to-use' sparse matrix format
		return mat.to_csr();
	}

	std::vector<double> approximate_rhs2() const {
		std::vector<double> ret(grid.n_points());

		//unsigned long a = grid.point(grid.n_points() - 1).x() / _hx;
		//unsigned long b = grid.point(grid.n_points() - 1).y() / _hy;

		//ret[0] = exact_solution(grid.point(0).x(), grid.point(0).y());
		// boundary left
		for (size_t i = 0; i < _ny; ++i)
		{
			ret[grid.to_linear_point_index({ 0, i })] = exact_solution(grid.point(0).x(), grid.point(i).y());
		}
		// boundary right
		for (size_t i = 0; i < _ny; ++i)
		{
			ret[grid.to_linear_point_index({ _nx - 1, i})] = exact_solution(grid.point(_nx - 1).x(), grid.point(i).y());
		}
		// boundary top
		for (size_t i = 0; i < _nx; ++i)
		{
			ret[grid.to_linear_point_index({ i, 0 })] = exact_solution(grid.point(i).x(), grid.point(0).y());
		}
		// boundary bottom
		for (size_t i = 0; i < _nx; ++i)
		{
			ret[grid.to_linear_point_index({ i, _ny - 1})] = exact_solution(grid.point(i).x(), grid.point(_ny - 1).y());
		}
		
		for (size_t i = 1; i < _ny - 1; ++i) {
			for (size_t j = 1; j < _nx - 1; ++j)
			{
				ret[grid.to_linear_point_index({ j, i})] = exact_rhsX(grid.point(grid.to_linear_point_index({ j, i})).x(), grid.point(grid.to_linear_point_index({ j, i})).y()) + 
					exact_rhsY(grid.point(grid.to_linear_point_index({ j, i})).x(), grid.point(grid.to_linear_point_index({ j, i})).y());
			}
		}
		//ret[grid.n_points() - 1] = exact_solution(grid.point(grid.n_points() - 1).x(), grid.point(grid.n_points() - 1).y());
		return ret;
	}

	double compute_norm2D() const {
		//size_t a = grid.point(grid.n_points() - 1).x() / (grid.point(1).x() - grid.point(0).x());
		//size_t b = grid.point(grid.n_points() - 1).y() / (grid.point(2 * a + 1).y() - grid.point(a).y());

		// weights
		//double hx = grid.point(1).x() - grid.point(0).x();
		//double hy = grid.point(2 * a + 1).y() - grid.point(a).y();
		std::vector<double> w(grid.n_points(), _hx * _hy);
		//unsigned long long ny = grid.point(grid.n_points() - 1).y() / _hy;
		//unsigned long long nx = grid.point(grid.n_points() - 1).x() / _hx;
		w[grid.to_linear_point_index({0, 0})] = w[grid.to_linear_point_index({ _nx - 1, 0 })] = w[grid.to_linear_point_index({ _nx - 1,  _ny - 1 })] = 
			w[grid.to_linear_point_index({ 0,  _ny - 1 })] = (_hx * _hy) / 4.0;
		for (size_t i = 1; i < _nx - 1; ++i)
		{
			w[grid.to_linear_point_index({ i, 0 })] = (_hx * _hy) / 2.0;
		}
		for (size_t i = 1; i < _ny - 1; ++i)
		{
			w[grid.to_linear_point_index({ _nx - 1, i })] = (_hx * _hy) / 2.0;
		}
		for (size_t i = 1; i < _nx - 1; ++i)
		{
			w[grid.to_linear_point_index({ i, _ny - 1 })] = (_hx * _hy) / 2.0;
		}
		for (size_t i = 1; i < _ny - 1; ++i)
		{
			w[grid.to_linear_point_index({ 0, i })] = (_hx * _hy) / 2.0;
		}

		// sum
		double sum = 0;
		for (size_t i = 0; i < grid.n_points(); ++i) {
			double diff = u[i] - exact_solution(grid.point(i).x(), grid.point(i).y());
			sum += w[i] * diff * diff;
		}

		double lenX = grid.point(grid.n_points() - 1).x() - grid.point(0).x();
		double lenY = grid.point(grid.n_points() - 1).y() - grid.point(0).y();
		return std::sqrt(sum / (lenX * lenY));
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

	// precalculated norm2 results for some n_cells values
	// used for CHECK procedures
	std::map<size_t, double> norm2_for_compare{
		/*{20, 0.179124},
		{200, 0.00158055},
		{2000, 1.57849e-05},*/
	};

	// loop over n_cells value
	//for (size_t n_cells : {10, 20, 50, 100, 200, 500, 1000}) {
	size_t n_cells = 3;
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
			auto found = norm2_for_compare.find(n_cells);
			if (found != norm2_for_compare.end()) {
				CHECK(n2 == Approx(found->second).margin(1e-6));
			}
	//}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////