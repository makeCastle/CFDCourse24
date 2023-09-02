#ifndef CFD24_GEOM_PRIMITIVES_HPP
#define CFD24_GEOM_PRIMITIVES_HPP

#include <array>

namespace cfd{

class Point: public std::array<double, 3>{
public:
	Point(double x=0, double y=0, double z=0): std::array<double, 3>({x, y, z}){}
};

}
#endif
