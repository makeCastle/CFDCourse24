#include "cfd24_test.hpp"
#include "cfd24/grid/grid1d.hpp"
#include "cfd24/grid/vtk.hpp"


TEST_CASE("Grid1d", "[grid1]"){
	cfd::Grid1D grid(0, 1, 3);

	CHECK(grid.n_points() == 4);
	CHECK(grid.n_cells() == 3);
	CHECK(grid.n_faces() == 4);

	grid.save_vtk("out.vtk");
	cfd::VtkUtils::add_point_data(std::vector<double>{0, 1, 2, 3}, "data1", "out.vtk");
}
