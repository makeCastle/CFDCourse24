﻿#include "cfd24_test.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/grid/unstructured_grid2d.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/mat/csrmat.hpp"
#include "cfd24/mat/lodmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/debug/printer.hpp"
#include "test/utils/filesystem.hpp"

using namespace cfd;

struct TestPoisson2FvmWorker{
	static double exact_solution(Point p){
		double x = p.x();
		double y = p.y();
		return cos(10*x*x)*sin(10*y) + sin(10*x*x)*cos(10*x);
	}
	static double exact_rhs(Point p){
		double x = p.x();
		double y = p.y();
		return (20*sin(10*x*x)+(400*x*x+100)*cos(10*x*x))*sin(10*y)
				+(400*x*x+100)*cos(10*x)*sin(10*x*x)
				+(400*x*sin(10*x)-20*cos(10*x))*cos(10*x*x);
	}

	TestPoisson2FvmWorker(const IGrid& grid);
	double solve();
	void save_vtk(const std::string& filename) const;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// partial derivative on x
	static double exact_dudx(Point p) {
		double x = p.x();
		double y = p.y();

		return (-20.0 * x * sin(10.0 * y) * sin(10 * x*x) + 20.0 * x * cos(10.0 * x * x) * cos(10.0 * x) - 10.0 * sin(10.0 * x* x) * sin(10.0 * x));
	}
	// partial derivative on y
	static double exact_dudy(Point p) {
		double x = p.x();
		double y = p.y();

		return (10.0 * cos(10.0 * x*x) * cos(10.0 * y));
	}
	static double exact_dudn(Point p) {
		double x = p.x();
		double y = p.y();
		if (std::abs(x) < 1e-6) {
			// left boundary -du/dx
			return -exact_dudx(p);
		}
		else if (std::abs(x - 1) < 1e-6) {
			// right boundary du/dx
			return exact_dudx(p);
		}
		else if (std::abs(y) < 1e-6) {
			// bottom boundary -du/dy
			return -exact_dudy(p);
		}
		else if (std::abs(y - 1) < 1e-6) {
			// top boundary du/dy
			return exact_dudy(p);
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

private:
	const IGrid& _grid;
	std::vector<size_t> _internal_faces;
	struct DirichletFace{
		size_t iface;
		size_t icell;
		double value;
		Vector outer_normal;
	};
	std::vector<DirichletFace> _dirichlet_faces;
	std::vector<double> _u;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	struct NeumannFace {
		size_t iface;
		size_t icell;
		double value;
	};
	std::vector<NeumannFace> _neumann_faces;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////


	CsrMatrix approximate_lhs() const;
	std::vector<double> approximate_rhs() const;
	double compute_norm2() const;
};

TestPoisson2FvmWorker::TestPoisson2FvmWorker(const IGrid& grid): _grid(grid){
	// assemble face lists
	for (size_t iface=0; iface<_grid.n_faces(); ++iface){
		size_t icell_negative = _grid.tab_face_cell(iface)[0];
		size_t icell_positive = _grid.tab_face_cell(iface)[1];
		if (icell_positive != INVALID_INDEX && icell_negative != INVALID_INDEX){
			// internal faces list
			_internal_faces.push_back(iface);
		} else {
			//// dirichlet faces list
			//DirichletFace dir_face;
			//dir_face.iface = iface;
			//dir_face.value = exact_solution(_grid.face_center(iface));
			//if (icell_positive == INVALID_INDEX){
			//	dir_face.icell = icell_negative;
			//	dir_face.outer_normal = _grid.face_normal(iface);
			//} else {
			//	dir_face.icell = icell_positive;
			//	dir_face.outer_normal = -_grid.face_normal(iface);
			//}
			
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//neumann faces list
			NeumannFace neum_face;
			neum_face.iface = iface;
			neum_face.value = exact_dudn(_grid.face_center(iface));
			if (icell_positive == INVALID_INDEX){
				neum_face.icell = icell_negative;
			} else {
				neum_face.icell = icell_positive;
			}
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//_dirichlet_faces.push_back(dir_face);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			_neumann_faces.push_back(neum_face);
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}
}

double TestPoisson2FvmWorker::solve(){
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
void TestPoisson2FvmWorker::save_vtk(const std::string& filename) const{
	// save grid
	_grid.save_vtk(filename);
	
	// save numerical solution
	VtkUtils::add_cell_data(_u, "numerical", filename);

	// save exact solution
	std::vector<double> exact(_grid.n_cells());
	for (size_t i=0; i<_grid.n_cells(); ++i){
		exact[i] = exact_solution(_grid.cell_center(i));
	}
	VtkUtils::add_cell_data(exact, "exact", filename);
}

CsrMatrix TestPoisson2FvmWorker::approximate_lhs() const{
	LodMatrix mat(_grid.n_cells());
	// internal faces
	for (size_t iface: _internal_faces){
		Vector normal = _grid.face_normal(iface);
		size_t negative_side_cell = _grid.tab_face_cell(iface)[0];
		size_t positive_side_cell = _grid.tab_face_cell(iface)[1];
		Point ci = _grid.cell_center(negative_side_cell);
		Point cj = _grid.cell_center(positive_side_cell);
		double h = dot_product(normal, cj-ci);
		double coef = _grid.face_area(iface) / h;

		mat.add_value(negative_side_cell, negative_side_cell, coef);
		mat.add_value(positive_side_cell, positive_side_cell, coef);
		mat.add_value(negative_side_cell, positive_side_cell, -coef);
		mat.add_value(positive_side_cell, negative_side_cell, -coef);
	}
	// dirichlet faces
	/*for (const DirichletFace& dir_face: _dirichlet_faces){
		size_t icell = dir_face.icell;
		size_t iface = dir_face.iface;
		Point gs = _grid.face_center(iface);
		Point ci = _grid.cell_center(icell);
		Vector normal = dir_face.outer_normal;
		double h = dot_product(normal, gs-ci);
		double coef = _grid.face_area(iface) / h;
		mat.add_value(icell, icell, coef);
	}*/

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	mat.set_unit_row(0);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return mat.to_csr();
}

std::vector<double> TestPoisson2FvmWorker::approximate_rhs() const{
	std::vector<double> rhs(_grid.n_cells(), 0.0);
	// internal
	for (size_t icell=0; icell < _grid.n_cells(); ++icell){
		double value = exact_rhs(_grid.cell_center(icell));
		double volume = _grid.cell_volume(icell);
		rhs[icell] = value * volume;
	}
	// dirichlet faces
	/*for (const DirichletFace& dir_face: _dirichlet_faces){
		size_t icell = dir_face.icell;
		size_t iface = dir_face.iface;
		Point gs = _grid.face_center(iface);
		Point ci = _grid.cell_center(icell);
		Vector normal = dir_face.outer_normal;
		double h = dot_product(normal, gs-ci);
		double coef = _grid.face_area(iface) / h;
		rhs[icell] += dir_face.value * coef;
	}*/

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// neumann faces
	for (const NeumannFace& neum_face : _neumann_faces) {
		size_t icell = neum_face.icell;
		size_t iface = neum_face.iface;
		rhs[icell] += _grid.face_area(iface) * neum_face.value; // neum_face.value: q = du/dn
	}

	rhs[0] = exact_solution(_grid.cell_center(0));
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return rhs;
}

double TestPoisson2FvmWorker::compute_norm2() const{
	double norm2 = 0;
	double full_area = 0;
	for (size_t icell=0; icell<_grid.n_cells(); ++icell){
		double diff = _u[icell] - exact_solution(_grid.cell_center(icell));
		norm2 += _grid.cell_volume(icell) * diff * diff;
		full_area += _grid.cell_volume(icell);
	}
	return std::sqrt(norm2/full_area);
}

TEST_CASE("Poisson-fvm 2D solver", "[poisson2-fvm]"){
	std::cout << std::endl << "--- cfd24_test [poisson2-fvm] --- " << std::endl;

	size_t nx = 150;
	RegularGrid2D grid(0.0, 1.0, 0.0, 1.0, nx, nx);
	TestPoisson2FvmWorker worker(grid);
	double nrm = worker.solve();
	worker.save_vtk("poisson2_fvm.vtk");
	std::cout << grid.n_cells() << " " << nrm << std::endl;

	CHECK(nrm == Approx(0.04371).margin(1e-4));
}

TEST_CASE("Poisson-fvm 2D solver pebigrid", "[poisson2-fvm-pebigrid]") {
	std::cout << std::endl << "--- cfd24_test [poisson2-fvm-pebigrid] --- " << std::endl;

	std::string fn = test_directory_file("pebigrid.vtk");
	UnstructuredGrid2D grid = UnstructuredGrid2D::vtk_read(fn);
	TestPoisson2FvmWorker worker(grid);

	//size_t nx = 20;
	//RegularGrid2D grid(0.0, 1.0, 0.0, 1.0, nx, nx);
	//TestPoisson2FvmWorker worker(grid);
	double nrm = worker.solve();
	worker.save_vtk("poisson2_fvm_pebigrid1.vtk");
	std::cout << grid.n_cells() << " " << nrm << std::endl;

	CHECK(nrm == Approx(0.04371).margin(1e-4));
}

TEST_CASE("Poisson-fvm 2D solver tetragrid", "[poisson2-fvm-tetragrid]") {
	std::cout << std::endl << "--- cfd24_test [poisson2-fvm-tetragrid] --- " << std::endl;

	std::string fn = test_directory_file("tetragrid.vtk");
	UnstructuredGrid2D grid = UnstructuredGrid2D::vtk_read(fn);
	TestPoisson2FvmWorker worker(grid);

	//size_t nx = 20;
	//RegularGrid2D grid(0.0, 1.0, 0.0, 1.0, nx, nx);
	//TestPoisson2FvmWorker worker(grid);
	double nrm = worker.solve();
	worker.save_vtk("poisson2_fvm_tetragrid.vtk");
	std::cout << grid.n_cells() << " " << nrm << std::endl;

	//CHECK(nrm == Approx(0.04371).margin(1e-4));
}