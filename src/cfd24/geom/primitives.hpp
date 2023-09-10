#ifndef CFD24_GEOM_PRIMITIVES_HPP
#define CFD24_GEOM_PRIMITIVES_HPP

#include <array>

namespace cfd{

/**
 * @brief 3D point
 *
 * Access to coordinate values can be achieved via operator[] or x/y/z() functions
 */
class Point: public std::array<double, 3>{
public:
	Point(double x=0, double y=0, double z=0): std::array<double, 3>({x, y, z}){}

	double x() const { return operator[](0); }
	double y() const { return operator[](1); }
	double z() const { return operator[](2); }

	double& x(){ return operator[](0); }
	double& y(){ return operator[](1); }
	double& z(){ return operator[](2); }
};

}
#endif
