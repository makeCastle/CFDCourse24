#include "cfd24_test.hpp"
#include "test/utils/filesystem.hpp"
#include "cfd24/grid/unstructured_grid2d.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/grid/grid1d.hpp"
#include "cfd24/mat/csrmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/fem/fem_assembler.hpp"
#include "cfd24/fem/fem_sorted_cell_info.hpp"
#include "cfd24/fem/elem1d/segment_linear.hpp"
#include "cfd24/fem/elem1d/segment_cubic.hpp"
#include "cfd24/fem/elem2d/triangle_linear.hpp"
#include "cfd24/fem/elem2d/triangle_cubic.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/debug/printer.hpp"
#include "cfd24/debug/saver.hpp"
#include "cfd24/fem/fem_numeric_integrals.hpp"
#include "cfd24/numeric_integration/segment_quadrature.hpp"
#include "cfd24/numeric_integration/triangle_quadrature.hpp"

using namespace cfd;

namespace {

///////////////////////////////////////////////////////////////////////////////
// Func
///////////////////////////////////////////////////////////////////////////////
struct Func1D: public IPointFunction{
	double value(Point p) const override{
		double x = p.x();
		return sin(10*x*x);
	}
	Vector grad(Point p) const override{
		double x = p.x();
		return { 20*x*cos(10*x*x) };
	};
} func1d;

struct Rhs1D: public IPointFunction{
	double value(Point p) const override{
		double x = p.x();
		return 400*x*x*sin(10*x*x) - 20*cos(10*x*x);
	};
	Vector grad(Point p) const override{
		double x = p.x();
		return { 1200*x*sin(10*x*x)+8000*x*x*x*cos(10*x*x) };
	};
} rhs1d;

struct Func2D: public IPointFunction{
	double value(Point p) const override{
		double x = p.x();
		double y = p.y();
		return cos(10*x*x)*sin(10*y) + sin(10*x*x)*cos(10*x);
	}
	Vector grad(Point p) const override{
		double x = p.x();
		double y = p.y();
		return {
			-(20*x*sin(10*x*x)*sin(10*y))-10*sin(10*x)*sin(10*x*x)+20*x*cos(10*x)*cos(10*x*x),
			10*cos(10*x*x)*cos(10*y)
		};
	};
} func2d;

struct Rhs2D: public IPointFunction{
	double value(Point p) const override{
		double x = p.x();
		double y = p.y();
		return (20*sin(10*x*x)+(400*x*x+100)*cos(10*x*x))*sin(10*y)
				+(400*x*x+100)*cos(10*x)*sin(10*x*x)
				+(400*x*sin(10*x)-20*cos(10*x))*cos(10*x*x);
	}
	Vector grad(Point p) const override{
		double x = p.x();
		double y = p.y();
		return {
			((-(8000*x*x*x)-2000*x)*sin(10*x*x)+1200*x*cos(10*x*x))*sin(10*y)
				+((-(12000*x*x)-1000)*sin(10*x)+1200*x*cos(10*x))*sin(10*x*x)
				+(600*sin(10*x)+(8000*x*x*x+6000*x)*cos(10*x))*cos(10*x*x),
			(200*sin(10*x*x)+(4000*x*x+1000)*cos(10*x*x))*cos(10*y)
		};
	};
} rhs2d;
}

///////////////////////////////////////////////////////////////////////////////
// ITestPoissonHermiteWorker
///////////////////////////////////////////////////////////////////////////////
struct ITestPoissonHermiteWorker{
	ITestPoissonHermiteWorker(const IGrid& grid, const IPointFunction& solution_func, const IPointFunction& rhs_func, const FemAssembler& fem);
	virtual double solve();
	void save_vtk(const std::string& filename) const;
protected:
	const IGrid& _grid;
	const IPointFunction& _solution_func;
	const IPointFunction& _rhs_func;
	FemAssembler _fem;

	std::vector<double> _u;

	CsrMatrix approximate_lhs() const;
	std::vector<double> approximate_rhs() const;
	double compute_norm2() const;
};

ITestPoissonHermiteWorker::ITestPoissonHermiteWorker(const IGrid& grid,
		const IPointFunction& solution_func,
		const IPointFunction& rhs_func,
		const FemAssembler& fem): _grid(grid), _solution_func(solution_func), _rhs_func(rhs_func), _fem(fem){
}

double ITestPoissonHermiteWorker::solve(){
	// 1. build SLAE
	CsrMatrix mat = approximate_lhs();
	std::vector<double> rhs = approximate_rhs();
	// 2. Dirichlet bc
	for (size_t ibas: _grid.boundary_points()){
		mat.set_unit_row(ibas);
		Point p = _fem.reference_point(ibas);
		rhs[ibas] = _solution_func(p);
	}
	// 3. solve SLAE
	AmgcMatrixSolver solver({ {"precond.relax.type", "gauss_seidel"} });
	solver.set_matrix(mat);
	solver.solve(rhs, _u);
	// 4. compute norm2
	return compute_norm2();
}

void ITestPoissonHermiteWorker::save_vtk(const std::string& filename) const{
	// save grid
	_grid.save_vtk(filename);
	// save numerical solution
	VtkUtils::add_point_data(_u, "numerical", filename, _grid.n_points());
	// save exact solution
	std::vector<double> exact(_grid.n_points());
	for (size_t i=0; i<_grid.n_points(); ++i){
		exact[i] = _solution_func(_grid.point(i));
	}
	VtkUtils::add_point_data(exact, "exact", filename, _grid.n_points());
}

CsrMatrix ITestPoissonHermiteWorker::approximate_lhs() const{
	CsrMatrix ret(_fem.stencil());
	for (size_t ielem=0; ielem < _fem.n_elements(); ++ielem){
		const FemElement& elem = _fem.element(ielem);
		std::vector<double> local_stiff = elem.integrals->stiff_matrix();
		_fem.add_to_global_matrix(ielem, local_stiff, ret.vals());
	}
	return ret;
}

std::vector<double> ITestPoissonHermiteWorker::approximate_rhs() const{
	// mass matrix
	CsrMatrix mass(_fem.stencil());
	for (size_t ielem=0; ielem < _fem.n_elements(); ++ielem){
		const FemElement& elem = _fem.element(ielem);
		std::vector<double> local_mass = elem.integrals->mass_matrix();
		_fem.add_to_global_matrix(ielem, local_mass, mass.vals());
	}
	// rhs = Mass * f
	std::vector<double> fvec = _fem.approximate(_rhs_func);
	return mass.mult_vec(fvec);
}

double ITestPoissonHermiteWorker::compute_norm2() const{
	// mass matrix
	CsrMatrix mass(_fem.stencil());
	for (size_t ielem=0; ielem < _fem.n_elements(); ++ielem){
		const FemElement& elem = _fem.element(ielem);
		std::vector<double> local_mass = elem.integrals->mass_matrix();
		_fem.add_to_global_matrix(ielem, local_mass, mass.vals());
	}
	std::vector<double> diff;
	std::vector<double> exact_u = _fem.approximate(_solution_func);
	for (size_t ibas=0; ibas<_fem.n_bases(); ++ibas){
		diff.push_back(_u[ibas] - exact_u[ibas]);
	}
	std::vector<double> diff_sum1 = mass.mult_vec(diff);
	double integral = 0;
	for (size_t i=0; i<diff_sum1.size(); ++i){
		integral += diff_sum1[i] * diff[i];
	}
	double full_area = 1;
	return std::sqrt(integral/full_area);
}

///////////////////////////////////////////////////////////////////////////////
// TestPoissonHermiteSegmentWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonHermiteSegmentWorker: public ITestPoissonHermiteWorker{
	TestPoissonHermiteSegmentWorker(const IGrid& grid): ITestPoissonHermiteWorker(grid, func1d, rhs1d, build_fem(grid)){ }
	static FemAssembler build_fem(const IGrid& grid);
};

FemAssembler TestPoissonHermiteSegmentWorker::build_fem(const IGrid& grid){
	std::vector<FemElement> elements;
	std::vector<std::vector<size_t>> tab_elem_basis;

	// elements
	for (size_t icell=0; icell < grid.n_cells(); ++icell){
		std::vector<size_t> ipoints = grid.tab_cell_point(icell);
		Point p0 = grid.point(ipoints[0]);
		Point p1 = grid.point(ipoints[1]);
		
		auto geom = std::make_shared<SegmentLinearGeometry>(p0, p1);
		auto basis = std::make_shared<SegmentHermiteBasis>(geom);
		auto integrals = std::make_shared<NumericElementIntegrals>(quadrature_segment_gauss4(), geom, basis);
		FemElement elem{geom, basis, integrals};

		elements.push_back(elem);
		tab_elem_basis.push_back({
			ipoints[0], ipoints[1],
			grid.n_points() + ipoints[0], grid.n_points() + ipoints[1]
		});
	}

	size_t n_bases = 2*grid.n_points();
	return FemAssembler(n_bases, elements, tab_elem_basis);
}

TEST_CASE("Poisson-fem 1D solver, hermite segments", "[poisson1-fem-hermite]"){
	std::cout << std::endl << "--- cfd24_test [poisson1-fem-hermite] --- " << std::endl;
	Grid1D grid(0, 1, 10);
	TestPoissonHermiteSegmentWorker worker(grid);
	double nrm = worker.solve();
	std::cout << grid.n_cells() << " " << nrm << std::endl;
	CHECK(nrm == Approx(0.00500192).margin(1e-6));

};

///////////////////////////////////////////////////////////////////////////////
// TestPoissonCubicHermiteWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonHermiteTriangleWorker: public ITestPoissonHermiteWorker{
	static FemAssembler build_fem(const IGrid& grid);
	TestPoissonHermiteTriangleWorker(const IGrid& grid): ITestPoissonHermiteWorker(grid, func2d, rhs2d, build_fem(grid)){ }

};

FemAssembler TestPoissonHermiteTriangleWorker::build_fem(const IGrid& grid){
	std::vector<FemElement> elements;
	std::vector<std::vector<size_t>> tab_elem_basis;

	// elements
	for (size_t icell=0; icell < grid.n_cells(); ++icell){
		PolygonElementInfo cell_info(grid, icell);

		if (cell_info.n_points() == 3){
			Point p0 = grid.point(cell_info.ipoints[0]);
			Point p1 = grid.point(cell_info.ipoints[1]);
			Point p2 = grid.point(cell_info.ipoints[2]);
			
			auto geom = std::make_shared<TriangleLinearGeometry>(p0, p1, p2);
			auto basis = std::make_shared<TriangleHermiteBasis>(geom);
			auto integrals = std::make_shared<NumericElementIntegrals>(quadrature_triangle_gauss6(), geom, basis);

			elements.push_back(FemElement{geom, basis, integrals});

			size_t bas3 = grid.n_points() + cell_info.ipoints[0];
			size_t bas4 = grid.n_points() + cell_info.ipoints[1];
			size_t bas5 = grid.n_points() + cell_info.ipoints[2];
			size_t bas6 = 2*grid.n_points() + cell_info.ipoints[0];
			size_t bas7 = 2*grid.n_points() + cell_info.ipoints[1];
			size_t bas8 = 2*grid.n_points() + cell_info.ipoints[2];
			size_t bas9 = 3*grid.n_points() + icell;

			tab_elem_basis.push_back({
				cell_info.ipoints[0], cell_info.ipoints[1], cell_info.ipoints[2],
				bas3, bas4, bas5, bas6, bas7, bas8,
				bas9
			});
		} else {
			throw std::runtime_error("Invalid fem grid");
		}
	}

	size_t n_bases = 3*grid.n_points() + grid.n_cells();

	return FemAssembler(n_bases, elements, tab_elem_basis);
}

TEST_CASE("Poisson-fem 2D solver, cubic hermite", "[poisson2-fem-hermite]"){
	std::cout << std::endl << "--- cfd24_test [poisson2-fem-hermite] --- " << std::endl;

	std::string grid_fn = test_directory_file("trigrid_500.vtk");
	//for (std::string s: {"500", "1k", "2k", "5k", "10k", "20k", "50k", "100k"}){
	//for (std::string s: {"500"}){
		//std::string grid_fn = "trigrid_" + s + ".vtk";
		//grid_fn = tmp_directory_file(grid_fn);
		UnstructuredGrid2D grid = UnstructuredGrid2D::vtk_read(grid_fn, true);
		TestPoissonHermiteTriangleWorker worker(grid);
		double nrm = worker.solve();
		std::cout << std::sqrt(grid.n_cells()) << " " << nrm << std::endl;
	//}
	//CHECK(nrm == Approx(0.0638327072).margin(1e-6));
}
