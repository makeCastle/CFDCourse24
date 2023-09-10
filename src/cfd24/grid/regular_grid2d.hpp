#ifndef CFD24_GEOM_REGULAR_GRID2_HPP
#define CFD24_GEOM_REGULAR_GRID2_HPP

#include "cfd24/grid/i_grid.hpp"

namespace cfd{

/**
 * @brief structured 2D quadraliteral grid with
 *
 * Fast indexing for points and cells goes along x dimension, so that the linear index equals
 * \f$
 * 	i = i_x + n_x i_y
 * \f$
 */
class RegularGrid2D: public IGrid2D{
public:
	/**
	 * @brief Split index as the [ix, iy] pair
	 */
	using split_index_t = std::array<size_t, 2>;

	/**
	 * @brief Builds regular equidistant 1d grid
	 *
	 * @param x0  x minimum
	 * @param x1  x maximum
	 * @param y0  y minimum
	 * @param y1  y maximum
	 * @param nx  number of cells in x direction
	 * @param ny  number of cells in y direction
	 */
	RegularGrid2D(double x0, double x1, double y0, double y1, size_t nx, size_t ny);

	/**
	 * @brief Converts split [ix, iy] point index to linear point index
	 *
	 * @param point_split_index  split point index as [ix, iy] array
	 * @return linear point index
	 */
	size_t to_linear_point_index(split_index_t point_split_index) const;

	/**
	 * @brief Converts linear point index to [ix, iy] point index
	 *
	 * @param ipoint  linear point index in [0, n_points()-1]
	 * @return split point index as the [ix, iy] array
	 */
	split_index_t to_split_point_index(size_t ipoint) const;

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
