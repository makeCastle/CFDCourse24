#include "cfd24_test.hpp"
#include "cfd24/grid/grid1d.hpp"
#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/grid/vtk.hpp"

using namespace cfd;


TEST_CASE("Grid1d", "[grid1]"){
	Grid1D grid(0, 1, 3);

	CHECK(grid.n_points() == 4);
	CHECK(grid.n_cells() == 3);
	CHECK(grid.n_faces() == 4);

	grid.save_vtk("out1.vtk");
	VtkUtils::add_point_data(std::vector<double>{0, 1, 2, 3}, "data1", "out1.vtk");
}

TEST_CASE("RegularGrid2d", "[reggrid2]"){
	RegularGrid2D grid(0, 1, 1, 3, 3, 2);

	CHECK(grid.n_points() == 12);
	CHECK(grid.n_cells() == 6);

	grid.save_vtk("out2.vtk");
	VtkUtils::add_point_data(std::vector<double>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, "data1", "out2.vtk");
}
