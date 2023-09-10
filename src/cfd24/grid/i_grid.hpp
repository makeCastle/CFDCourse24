#ifndef CFD24_GRID_I_GRID_HPP
#define CFD24_GRID_I_GRID_HPP

#include "cfd24/cfd_common.hpp"
#include "cfd24/geom/primitives.hpp"

namespace cfd{

/**
 * @brief Any dimension grid interface
 */
class IGrid{
public:
	virtual ~IGrid() = default;

	/// @brief number of geometric dimensions: 1, 2 or 3
	virtual size_t dim() const = 0;

	/// @brief number of grid points
	virtual size_t n_points() const = 0;

	/// @brief number of grid cells
	virtual size_t n_cells() const = 0;

	/// @brief number of grid faces
	virtual size_t n_faces() const = 0;

	/// @brief get point at i-th index
	virtual Point point(size_t ipoint) const = 0;

	/// @brief get list of all points
	virtual std::vector<Point> points() const = 0;

	/**
	 * @brief cell->point connectivity table
	 *
	 * @param icell  cell index
	 *
	 * For 1d grid points are ordered by x coordinate \n
	 * For 2d grid points have counter clockwise direction
	 */
	virtual std::vector<size_t> tab_cell_point(size_t icell) const = 0;

	/**
	 * @brief Saves grid to vtk format
	 *
	 * @param fname  output filename
	 */
	virtual void save_vtk(std::string fname) const = 0;
};


/**
 * @brief Abstract 1D grid
 */
class IGrid1D: public IGrid{
public:
	virtual ~IGrid1D() = default;

	size_t dim() const override { return 1; }
};

/**
 * @brief Abstract 2D grid
 */
class IGrid2D: public IGrid{
public:
	virtual ~IGrid2D() = default;

	size_t dim() const override { return 2; }
};

/**
 * @brief Abstract 3D grid
 */
class IGrid3D: public IGrid{
public:
	virtual ~IGrid3D() = default;

	size_t dim() const override { return 3; }
};


}
#endif
