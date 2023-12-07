#ifndef CFD24_GEOM_POINT_HPP
#define CFD24_GEOM_POINT_HPP

#include <array>
#include <cmath>

namespace cfd{

/**
 * @brief 3D point
 *
 * Access to coordinate values can be achieved via operator[] or x/y/z() methods
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

	Point& operator+=(const Point& other){
		x() += other.x();
		y() += other.y();
		z() += other.z();
		return *this;
	}
	friend Point operator+(const Point& p1, const Point& p2){
		Point p(p1);
		p += p2;
		return p;
	}

	Point& operator-=(const Point& other){
		x() -= other.x();
		y() -= other.y();
		z() -= other.z();
		return *this;
	}
	friend Point operator-(const Point& p1, const Point& p2){
		Point p(p1);
		p -= p2;
		return p;
	}
	friend Point operator-(const Point& p1){
		Point p(p1);
		p.x() *= -1;
		p.y() *= -1;
		p.z() *= -1;
		return p;
	}

	Point& operator*=(double a){
		x() *= a;
		y() *= a;
		z() *= a;
		return *this;
	}
	friend Point operator*(double a, const Point& p1){
		Point p(p1);
		p *= a;
		return p;
	}

	Point& operator/=(double a){
		x() /= a;
		y() /= a;
		z() /= a;
		return *this;
	}
	friend Point operator/(const Point& p1, double a){
		Point p(p1);
		p /= a;
		return p;
	}

};


/**
 * @brief 3D vector
 *
 * Alias to cfd::Point.
 * Access to coordinate values can be achieved via operator[] or x/y/z() methods
 */
using Vector = Point;

/**
 * @brief dot product of two vectors
 */
inline double dot_product(const Vector& v1, const Vector& v2){
	return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

/**
 * @brief cross product of two vectors;
 */
inline Vector cross_product(const Vector& v1, const Vector& v2){
	return Point(
		v1.y()*v2.z() - v1.z()*v2.y(),
		v1.z()*v2.x() - v1.x()*v2.z(),
		v1.x()*v2.y() - v1.y()*v2.x());
}

/**
 * @brief z component of the vector cross product
 */
inline double cross_product_2d(const Vector& v1, const Vector& v2){
	return v1.x()*v2.y() - v1.y()*v2.x();
}

/**
 * @brief vector length
 */
inline double vector_abs(const Vector& v){
	return std::sqrt(v.x()*v.x() + v.y()*v.y() + v.z()*v.z());
}

/**
 * @brief vector squared length
 */
inline double vector_meas(const Vector& v){
	return v.x()*v.x() + v.y()*v.y() + v.z()*v.z();
}

}
#endif
