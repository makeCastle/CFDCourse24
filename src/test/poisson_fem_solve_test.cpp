#include "cfd24_test.hpp"
#include "test/utils/filesystem.hpp"
#include "cfd24/grid/unstructured_grid2d.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/grid/grid1d.hpp"
#include "cfd24/mat/csrmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/fem/fem_assembler.hpp"
#include "cfd24/fem/elem1d/segment_linear.hpp"
#include "cfd24/fem/elem2d/triangle_linear.hpp"
#include "cfd24/fem/elem2d/quadrangle_linear.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/debug/printer.hpp"
#include "cfd24/fem/fem_numeric_integrals.hpp"
#include "cfd24/numeric_integration/segment_quadrature.hpp"
#include "cfd24/numeric_integration/square_quadrature.hpp"
#include "cfd24/numeric_integration/triangle_quadrature.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// ITestPoissonFemWorker
///////////////////////////////////////////////////////////////////////////////
struct ITestPoissonFemWorker {
	virtual double exact_solution(Point p) const = 0;
	virtual double exact_rhs(Point p) const = 0;
	virtual std::vector<size_t> dirichlet_bases() const = 0;

	ITestPoissonFemWorker(const IGrid& grid, const FemAssembler& fem);
	double solve();
	void save_vtk(const std::string& filename) const;
protected:
	const IGrid& _grid;
	FemAssembler _fem;
	std::vector<double> _u;

	CsrMatrix approximate_lhs() const;
	std::vector<double> approximate_rhs() const;
	double compute_norm2() const;
};

ITestPoissonFemWorker::ITestPoissonFemWorker(const IGrid& grid, const FemAssembler& fem) : _grid(grid), _fem(fem) {
}

double ITestPoissonFemWorker::solve() {
	// 1. build SLAE
	CsrMatrix mat = approximate_lhs();
	std::vector<double> rhs = approximate_rhs();
	// 2. Dirichlet bc
	for (size_t ibas : dirichlet_bases()) {
		mat.set_unit_row(ibas);
		Point p = _fem.reference_point(ibas);
		rhs[ibas] = exact_solution(p);
	}
	// 3. solve SLAE
	AmgcMatrixSolver solver;
	solver.set_matrix(mat);
	solver.solve(rhs, _u);
	// 4. compute norm2
	return compute_norm2();
}

void ITestPoissonFemWorker::save_vtk(const std::string& filename) const {
	// save grid
	_grid.save_vtk(filename);
	// save numerical solution
	VtkUtils::add_point_data(_u, "numerical", filename);
	// save exact solution
	std::vector<double> exact(_grid.n_points());
	for (size_t i = 0; i < _grid.n_points(); ++i) {
		exact[i] = exact_solution(_grid.point(i));
	}
	VtkUtils::add_point_data(exact, "exact", filename);
}

CsrMatrix ITestPoissonFemWorker::approximate_lhs() const {
	CsrMatrix ret(_fem.stencil());
	for (size_t ielem = 0; ielem < _fem.n_elements(); ++ielem) {
		const FemElement& elem = _fem.element(ielem);
		std::vector<double> local_stiff = elem.integrals->stiff_matrix();
		_fem.add_to_global_matrix(ielem, local_stiff, ret.vals());
	}
	return ret;
}

std::vector<double> ITestPoissonFemWorker::approximate_rhs() const {
	// mass matrix
	CsrMatrix mass(_fem.stencil());
	for (size_t ielem = 0; ielem < _fem.n_elements(); ++ielem) {
		const FemElement& elem = _fem.element(ielem);
		std::vector<double> local_mass = elem.integrals->mass_matrix();
		_fem.add_to_global_matrix(ielem, local_mass, mass.vals());
	}
	// rhs = Mass * f
	std::vector<double> fvec(_fem.n_bases());
	for (size_t ibas = 0; ibas < _fem.n_bases(); ++ibas) {
		Point p = _fem.reference_point(ibas);
		fvec[ibas] = exact_rhs(p);
	}
	return mass.mult_vec(fvec);
}

double ITestPoissonFemWorker::compute_norm2() const {
	std::vector<double> force_vec(_fem.n_bases(), 0);
	for (size_t ielem = 0; ielem < _fem.n_elements(); ++ielem) {
		const FemElement& elem = _fem.element(ielem);
		std::vector<double> v = elem.integrals->load_vector();
		_fem.add_to_global_vector(ielem, v, force_vec);
	}
	double integral = 0;
	double full_area = 0;
	for (size_t ipoint = 0; ipoint < _grid.n_points(); ++ipoint) {
		double diff = _u[ipoint] - exact_solution(_grid.point(ipoint));
		integral += force_vec[ipoint] * (diff * diff);
		full_area += force_vec[ipoint];
	}
	return std::sqrt(integral / full_area);
}

////////////////////////////////////////////////////////////////////////////////
// ITestPoisson1FemWorker
////////////////////////////////////////////////////////////////////////////////
struct ITestPoisson1FemWorker : public ITestPoissonFemWorker {
	double exact_solution(Point p) const override {
		double x = p.x();
		return sin(10 * x * x);
	}
	double exact_rhs(Point p) const override {
		double x = p.x();
		return 400 * x * x * sin(10 * x * x) - 20 * cos(10 * x * x);
	}
	ITestPoisson1FemWorker(const IGrid& grid, const FemAssembler& fem) : ITestPoissonFemWorker(grid, fem) { }
};

////////////////////////////////////////////////////////////////////////////////
// ITestPoisson2FemWorker
////////////////////////////////////////////////////////////////////////////////
struct ITestPoisson2FemWorker : public ITestPoissonFemWorker {
	double exact_solution(Point p) const override {
		double x = p.x();
		double y = p.y();
		return cos(10 * x * x) * sin(10 * y) + sin(10 * x * x) * cos(10 * x);
	}
	double exact_rhs(Point p) const override {
		double x = p.x();
		double y = p.y();
		return (20 * sin(10 * x * x) + (400 * x * x + 100) * cos(10 * x * x)) * sin(10 * y)
			+ (400 * x * x + 100) * cos(10 * x) * sin(10 * x * x)
			+ (400 * x * sin(10 * x) - 20 * cos(10 * x)) * cos(10 * x * x);
	}
	ITestPoisson2FemWorker(const IGrid& grid, const FemAssembler& fem) : ITestPoissonFemWorker(grid, fem) { }
};

///////////////////////////////////////////////////////////////////////////////
// TestPoissonLinearSegmentWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonLinearSegmentWorker : public ITestPoisson1FemWorker {
	std::vector<size_t> dirichlet_bases() const override {
		return _grid.boundary_points();
	}

	TestPoissonLinearSegmentWorker(const IGrid& grid) : ITestPoisson1FemWorker(grid, build_fem(grid)) { }

	static FemAssembler build_fem(const IGrid& grid);
};

FemAssembler TestPoissonLinearSegmentWorker::build_fem(const IGrid& grid) {
	size_t n_bases = grid.n_points();
	std::vector<FemElement> elements;
	std::vector<std::vector<size_t>> tab_elem_basis;

	// elements
	for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
		std::vector<size_t> ipoints = grid.tab_cell_point(icell);
		Point p0 = grid.point(ipoints[0]);
		Point p1 = grid.point(ipoints[1]);

		auto geom = std::make_shared<SegmentLinearGeometry>(p0, p1);
		auto basis = std::make_shared<SegmentLinearBasis>();
		auto integrals = std::make_shared<SegmentLinearIntegrals>(geom->jacobi({}));
		FemElement elem{ geom, basis, integrals };

		elements.push_back(elem);
		tab_elem_basis.push_back(ipoints);
	}

	return FemAssembler(n_bases, elements, tab_elem_basis);
}

TEST_CASE("Poisson-fem 1D solver, segments", "[poisson1-fem-segm]") {
	std::cout << std::endl << "--- cfd24_test [poisson1-fem-segm] --- " << std::endl;
	Grid1D grid(0, 1, 10);
	TestPoissonLinearSegmentWorker worker(grid);
	double nrm = worker.solve();
	worker.save_vtk("poisson1_fem.vtk");
	std::cout << grid.n_cells() << " " << nrm << std::endl;
	CHECK(nrm == Approx(0.138156).margin(1e-6));
}

///////////////////////////////////////////////////////////////////////////////
// TestPoissonLinearTriangleWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonLinearTriangleWorker : public ITestPoisson2FemWorker {
	std::vector<size_t> dirichlet_bases() const override {
		return _grid.boundary_points();
	}
	static FemAssembler build_fem(const IGrid& grid);
	TestPoissonLinearTriangleWorker(const IGrid& grid) : ITestPoisson2FemWorker(grid, build_fem(grid)) { }
};

FemAssembler TestPoissonLinearTriangleWorker::build_fem(const IGrid& grid) {
	size_t n_bases = grid.n_points();
	std::vector<FemElement> elements;
	std::vector<std::vector<size_t>> tab_elem_basis;

	// elements
	for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
		std::vector<size_t> ipoints = grid.tab_cell_point(icell);
		Point p0 = grid.point(ipoints[0]);
		Point p1 = grid.point(ipoints[1]);
		Point p2 = grid.point(ipoints[2]);

		auto geom = std::make_shared<TriangleLinearGeometry>(p0, p1, p2);
		auto basis = std::make_shared<TriangleLinearBasis>();
		JacobiMatrix jac = geom->jacobi({ 0, 0 });
		auto integrals = std::make_shared<TriangleLinearIntegrals>(jac);
		FemElement elem{ geom, basis, integrals };

		elements.push_back(elem);
		tab_elem_basis.push_back(ipoints);
	}

	return FemAssembler(n_bases, elements, tab_elem_basis);
}

TEST_CASE("Poisson-fem 2D solver, triangles", "[poisson2-fem-tri]") {
	std::cout << std::endl << "--- cfd24_test [poisson2-fem-tri] --- " << std::endl;

	std::string grid_fn = test_directory_file("trigrid.vtk");
	UnstructuredGrid2D grid = UnstructuredGrid2D::vtk_read(grid_fn);
	TestPoissonLinearTriangleWorker worker(grid);
	double nrm = worker.solve();
	worker.save_vtk("poisson2_fem.vtk");
	std::cout << grid.n_cells() << " " << nrm << std::endl;
	CHECK(nrm == Approx(0.0638327072).margin(1e-6));
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TestPoissonTriQuadrangleWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonLinearTriQuadWorker : public ITestPoisson2FemWorker {
	std::vector<size_t> dirichlet_bases() const override {
		return _grid.boundary_points();
	}
	static FemAssembler build_fem(const IGrid& grid);
	TestPoissonLinearTriQuadWorker(const IGrid& grid) :
		ITestPoisson2FemWorker(grid, build_fem(grid)) { }
};

FemAssembler TestPoissonLinearTriQuadWorker::build_fem(const IGrid& grid) {
	size_t n_bases = grid.n_points();
	std::vector<FemElement> elements;
	std::vector<std::vector<size_t>> tab_elem_basis;

	// elements
	for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
		std::vector<size_t> ipoints = grid.tab_cell_point(icell);

		if (ipoints.size() == 3)
		{
			Point p0 = grid.point(ipoints[0]);
			Point p1 = grid.point(ipoints[1]);
			Point p2 = grid.point(ipoints[2]);

			auto geom = std::make_shared<TriangleLinearGeometry>(p0, p1, p2);
			auto basis = std::make_shared<TriangleLinearBasis>();
			JacobiMatrix jac = geom->jacobi({ 0, 0 });
			auto integrals = std::make_shared<TriangleLinearIntegrals>(jac);
			FemElement elem{ geom, basis, integrals };
			elements.push_back(elem);
			tab_elem_basis.push_back(ipoints);
		}
		else {
			Point p0 = grid.point(ipoints[0]);
			Point p1 = grid.point(ipoints[1]);
			Point p2 = grid.point(ipoints[2]);
			Point p3 = grid.point(ipoints[3]);
			auto geom = std::make_shared<QuadrangleLinearGeometry>(p0, p1, p2, p3);
			auto basis = std::make_shared<QuadrangleLinearBasis>();
			const Quadrature* quadrature = quadrature_square_gauss2();
			auto integrals = std::make_shared<NumericElementIntegrals>(quadrature, geom, basis); //quadrature_square_gauss4()
			FemElement elem{ geom, basis, integrals };
			elements.push_back(elem);
			tab_elem_basis.push_back(ipoints);
		}
	}

	return FemAssembler(n_bases, elements, tab_elem_basis);
}

TEST_CASE("Poisson-fem 2D solver, triangles & quadrangles", "[poisson2-fem-tri-quad]") {
	std::cout << std::endl << "--- cfd24_test [poisson2-fem-tri-quad] --- " << std::endl;

	std::string grid_fn = test_directory_file("tetragrid.vtk");
	UnstructuredGrid2D grid = UnstructuredGrid2D::vtk_read(grid_fn);
	TestPoissonLinearTriQuadWorker worker(grid);
	double nrm = worker.solve();
	worker.save_vtk("poisson2_fem_tetra_quad.vtk");
	std::cout << grid.n_cells() << " " << nrm << std::endl;
	CHECK(nrm == Approx(0.0638327072).margin(1e-6));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// TestPoissonBilinearQuadrangleWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonBilinearQuadrangleWorker : public ITestPoisson2FemWorker {
	std::vector<size_t> dirichlet_bases() const override {
		return _grid.boundary_points();
	}
	static FemAssembler build_fem(const IGrid& grid);
	TestPoissonBilinearQuadrangleWorker(const IGrid& grid) : ITestPoisson2FemWorker(grid, build_fem(grid)) { }
};

FemAssembler TestPoissonBilinearQuadrangleWorker::build_fem(const IGrid& grid) {
	size_t n_bases = grid.n_points();
	std::vector<FemElement> elements;
	std::vector<std::vector<size_t>> tab_elem_basis;

	// elements
	for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
		std::vector<size_t> ipoints = grid.tab_cell_point(icell);
		Point p0 = grid.point(ipoints[0]);
		Point p1 = grid.point(ipoints[1]);
		Point p2 = grid.point(ipoints[2]);
		Point p3 = grid.point(ipoints[3]);

		auto geom = std::make_shared<QuadrangleLinearGeometry>(p0, p1, p2, p3);
		auto basis = std::make_shared<QuadrangleLinearBasis>();
		const Quadrature* quadrature = quadrature_square_gauss2();
		auto integrals = std::make_shared<NumericElementIntegrals>(quadrature, geom, basis);
		FemElement elem{ geom, basis, integrals };

		elements.push_back(FemElement{ geom, basis, integrals });
		tab_elem_basis.push_back(ipoints);
	}

	return FemAssembler(n_bases, elements, tab_elem_basis);
}

TEST_CASE("Poisson-fem 2D solver, quadrangles", "[poisson2-fem-quad]") {
	std::cout << std::endl << "--- cfd24_test [poisson2-fem-quad] --- " << std::endl;
	RegularGrid2D grid(0, 1, 0, 1, 10, 10);
	TestPoissonLinearTriangleWorker worker(grid);
	double nrm = worker.solve();
	worker.save_vtk("poisson2_fem.vtk");
	std::cout << grid.n_cells() << " " << nrm << std::endl;
	CHECK(nrm == Approx(0.1585788505).margin(1e-6));
}

///////////////////////////////////////////////////////////////////////////////
// TestPoissonLinearTriangleGaussWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonLinearTriangleGaussWorker : public ITestPoisson2FemWorker {
	std::vector<size_t> dirichlet_bases() const override {
		return _grid.boundary_points();
	}
	static FemAssembler build_fem(const IGrid& grid);
	TestPoissonLinearTriangleGaussWorker(const IGrid& grid) : ITestPoisson2FemWorker(grid, build_fem(grid)) { }
};

FemAssembler TestPoissonLinearTriangleGaussWorker::build_fem(const IGrid& grid) {
	size_t n_bases = grid.n_points();
	std::vector<FemElement> elements;
	std::vector<std::vector<size_t>> tab_elem_basis;

	// elements
	for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
		std::vector<size_t> ipoints = grid.tab_cell_point(icell);
		Point p0 = grid.point(ipoints[0]);
		Point p1 = grid.point(ipoints[1]);
		Point p2 = grid.point(ipoints[2]);

		auto geom = std::make_shared<TriangleLinearGeometry>(p0, p1, p2);
		auto basis = std::make_shared<TriangleLinearBasis>();
		auto integrals = std::make_shared<NumericElementIntegrals>(quadrature_triangle_gauss4(), geom, basis);

		elements.push_back(FemElement{ geom, basis, integrals });
		tab_elem_basis.push_back(ipoints);
	}

	return FemAssembler(n_bases, elements, tab_elem_basis);
}

TEST_CASE("Poisson-fem 2D solver, triangles numerical integration", "[poisson2-fem-tri-quadrature]") {
	std::cout << std::endl << "--- cfd24_test [poisson2-fem-tri-quadrature] --- " << std::endl;

	std::string grid_fn = test_directory_file("trigrid_500.vtk");
	UnstructuredGrid2D grid = UnstructuredGrid2D::vtk_read(grid_fn, true);
	TestPoissonLinearTriangleGaussWorker worker(grid);
	double nrm = worker.solve();
	CHECK(nrm == Approx(0.0638327072).margin(1e-6));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TestPoissonTriQuadrangleQuadratureWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonLinearTriQuadQuadratureWorker : public ITestPoisson2FemWorker {
	std::vector<size_t> dirichlet_bases() const override {
		return _grid.boundary_points();
	}
	static FemAssembler build_fem(const IGrid& grid);
	TestPoissonLinearTriQuadQuadratureWorker(const IGrid& grid) :
		ITestPoisson2FemWorker(grid, build_fem(grid)) { }
};

FemAssembler TestPoissonLinearTriQuadQuadratureWorker::build_fem(const IGrid& grid) {
	size_t n_bases = grid.n_points();
	std::vector<FemElement> elements;
	std::vector<std::vector<size_t>> tab_elem_basis;

	// elements
	for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
		std::vector<size_t> ipoints = grid.tab_cell_point(icell);

		if (ipoints.size() == 3)
		{
			Point p0 = grid.point(ipoints[0]);
			Point p1 = grid.point(ipoints[1]);
			Point p2 = grid.point(ipoints[2]);

			auto geom = std::make_shared<TriangleLinearGeometry>(p0, p1, p2);
			auto basis = std::make_shared<TriangleLinearBasis>();
			auto integrals = std::make_shared<NumericElementIntegrals>(quadrature_triangle_gauss2(), geom, basis);

			elements.push_back(FemElement{ geom, basis, integrals });
			tab_elem_basis.push_back(ipoints);
		}
		else {
			Point p0 = grid.point(ipoints[0]);
			Point p1 = grid.point(ipoints[1]);
			Point p2 = grid.point(ipoints[2]);
			Point p3 = grid.point(ipoints[3]);
			auto geom = std::make_shared<QuadrangleLinearGeometry>(p0, p1, p2, p3);
			auto basis = std::make_shared<QuadrangleLinearBasis>();
			auto integrals = std::make_shared<NumericElementIntegrals>(quadrature_square_gauss2(), geom, basis);
			FemElement elem{ geom, basis, integrals };
			elements.push_back(elem);
			tab_elem_basis.push_back(ipoints);
		}
	}

	return FemAssembler(n_bases, elements, tab_elem_basis);
}
TEST_CASE("Poisson-fem 2D solver, triangles and quadrangles numerical integration", "[poisson2-fem-tri-quad-quadrature]") {
	std::cout << std::endl << "--- cfd24_test [poisson2-fem-tri-quad-quadrature] --- " << std::endl;

	std::string grid_fn = test_directory_file("tetragrid.vtk");
	UnstructuredGrid2D grid = UnstructuredGrid2D::vtk_read(grid_fn, true);
	TestPoissonLinearTriQuadQuadratureWorker worker(grid);
	double nrm = worker.solve();
	worker.save_vtk("poisson2_fem-tri4-quad.vtk");
	std::cout << grid.n_cells() << " " << nrm << std::endl;
	CHECK(nrm == Approx(0.0638327072).margin(1e-6));
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////