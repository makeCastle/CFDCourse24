#include "cfd24_test.hpp"
#include "cfd24/grid/grid1d.hpp"
#include "cfd24/mat/csrmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/mat/lodmat.hpp"
#include "C:/git_repos/CFDCourse24/src/cfd24/grid/regular_grid2d.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Poisson 1D solver
///////////////////////////////////////////////////////////////////////////////

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

	TestPoisson2Worker(size_t n_cellsX, size_t n_cellsY) : grid(0.0, 1.0, 0.0, 1.0, n_cellsX, n_cellsY) {}

	// returns norm2(u - u_exact)
	double solve() {
		// 1. build SLAE
		CsrMatrix mat = approximate_lhs2();
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

private:
	const RegularGrid2D grid;
	std::vector<double> u;

	CsrMatrix approximate_lhs2() const {
		// constant h = x[1] - x[0]
		double hx = grid.point(1).x() - grid.point(0).x();
		double hy = grid.point(1).y() - grid.point(0).y();

		// fill using 'easy-to-construct' sparse matrix format
		LodMatrix mat(grid.n_points());
		mat.add_value(0, 0, 1);
		mat.add_value(grid.n_points() - 1, grid.n_points() - 1, 1);
		double diag = 2.0 / (hx * hx) + 2.0 / (hy * hy);
		double nondiag1 = -1.0 / (hy * hy);
		double nondiag2 = -1.0 / (hx * hx);
		for (size_t i = 1; i < grid.n_points() - 1; ++i) {
			mat.add_value(i, i - 1, nondiag1);
			mat.add_value(i, i + 1, nondiag1);
			mat.add_value(i, i, diag);
		}

		// return 'easy-to-use' sparse matrix format
		return mat.to_csr();
	}

	std::vector<double> approximate_rhs2() const {
		std::vector<double> ret(grid.n_points());
		ret[0] = exact_solution(grid.point(0).x(), grid.point(0).y());
		ret[grid.n_points() - 1] = exact_solution(grid.point(grid.n_points() - 1).x(), grid.point(grid.n_points() - 1).y());
		for (size_t i = 1; i < grid.n_points() - 1; ++i) {
			for (size_t j = 0; j < grid.n_points() - 1; ++j)
			{
				ret[grid.to_linear_point_index({ i, j })] = exact_rhsX(grid.point(grid.to_linear_point_index({ i, j })).x(), grid.point(grid.to_linear_point_index({ i, j })).y()) + 
					exact_rhsY(grid.point(grid.to_linear_point_index({ i, j })).x(), grid.point(grid.to_linear_point_index({ i, j })).y());
			}
			
		}
		return ret;
	}

	double compute_norm2D() const {
		// weights
		double hx = grid.point(1).x() - grid.point(0).x();
		double hy = grid.point(1).y() - grid.point(0).y();
		std::vector<double> w(grid.n_points(), hx * hy);
		unsigned long long ny = grid.point(grid.n_points() - 1).y() / hy;
		unsigned long long nx = grid.point(grid.n_points() - 1).x() / hx;
		w[grid.to_linear_point_index({0, 0})] = w[grid.to_linear_point_index({ 0,  ny})] = w[grid.to_linear_point_index({ nx,  ny })] = 
			w[grid.to_linear_point_index({ nx,  0 })] = (hx * hy) / 4.0;
		for (size_t i = 1; i < nx; ++i)
		{
			w[grid.to_linear_point_index({ i, 0 })] = (hx * hy) / 2.0;
		}
		for (size_t i = 1; i < ny; ++i)
		{
			w[grid.to_linear_point_index({ nx, i })] = (hx * hy) / 2.0;
		}
		for (size_t i = 1; i < nx; ++i)
		{
			w[grid.to_linear_point_index({ i, ny })] = (hx * hy) / 2.0;
		}
		for (size_t i = 1; i < ny; ++i)
		{
			w[grid.to_linear_point_index({ 0, i })] = (hx * hy) / 2.0;
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

TEST_CASE("Poisson 2D solver", "[poisson2]") {
	std::cout << std::endl << "--- cfd24_test [poisson2] --- " << std::endl;

	// precalculated norm2 results for some n_cells values
	// used for CHECK procedures
	std::map<size_t, double> norm2_for_compare{
		{20, 0.179124},
		{200, 0.00158055},
		{2000, 1.57849e-05},
	};

	// loop over n_cells value
	for (size_t n_cellsX : {10, 20, 50, 100, 200, 500, 1000}) {
		// build test solver
		for (size_t n_cellsY : {10, 20, 50, 100, 200, 500, 1000}) {
			// build test solver
			TestPoisson2Worker worker(n_cellsX, n_cellsY);

			// solve and find norm2
			double n2 = worker.solve();

			// save into poisson1_ncells={n_cells}.vtk
			worker.save_vtk("poisson2_ncells=" + std::to_string(n_cellsX + n_cellsY) + ".vtk");

			// print (N_CELLS, NORM2) table entry
			std::cout << n_cellsX + n_cellsY << " " << n2 << std::endl;

			// CHECK if result for this n_cells
			// presents in the norm2_for_compare dictionary
			auto found = norm2_for_compare.find(n_cellsX + n_cellsY);
			if (found != norm2_for_compare.end()) {
				CHECK(n2 == Approx(found->second).margin(1e-6));
			}
		}
	}
}
