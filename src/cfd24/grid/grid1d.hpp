#ifndef CFD24_GEOM_GRID1_HPP
#define CFD24_GEOM_GRID1_HPP

#include "cfd24/grid/i_grid.hpp"

namespace cfd{

/**
 * 1D grid with ordered points
 */
class Grid1D: public IGrid1D{
public:
	/**
	 * Builds regular equidistant 1d grid
	 *
	 * @param left     left boundary
	 * @param right    right boundary
	 * @param n_cells  number of cells
	 */
	Grid1D(double left, double right, size_t n_cells);

	// overridden methods
	size_t n_points() const override;
	size_t n_cells() const override;
	size_t n_faces() const override;
	Point point(size_t ipoint) const override;
	std::vector<Point> points() const override;
	std::vector<size_t> tab_cell_point(size_t icell) const override;
	void save_vtk(std::string fname) const override;
private:
	std::vector<Point> _points;
};

}
#endif