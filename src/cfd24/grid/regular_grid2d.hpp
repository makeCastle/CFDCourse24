#ifndef CFD24_GEOM_REGULAR_GRID2_HPP
#define CFD24_GEOM_REGULAR_GRID2_HPP

#include "cfd24/grid/i_grid.hpp"

namespace cfd{

/**
 * 2D grid with ordered points
 */
class RegularGrid2D: public IGrid2D{
public:
	/**
	 * Builds regular equidistant 1d grid
	 *
	 * @param x0  x minimum
	 * @param x1  x maximum
	 * @param y0  y minimum
	 * @param y1  y maximum
	 * @param nx  number of cells in x direction
	 * @param ny  number of cells in y direction
	 */
	RegularGrid2D(double x0, double x1, double y0, double y1, size_t nx, size_t ny);

	// overridden methods
	size_t n_points() const override;
	size_t n_cells() const override;
	size_t n_faces() const override;
	Point point(size_t ipoint) const override;
	std::vector<Point> points() const override;
	std::vector<size_t> tab_cell_point(size_t icell) const override;
	void save_vtk(std::string fname) const override;
private:
	std::vector<double> _x;
	std::vector<double> _y;
};

}
#endif
