#ifndef CFD_SIMPLEX_HPP
#define CFD_SIMPLEX_HPP

#include "cfd24/geom/primitives.hpp"

namespace cfd{

/**
 * @brief signed triangle area
 */
double triangle_area(Point p0, Point p1, Point p2); 

}

#endif
