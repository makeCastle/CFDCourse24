#ifndef __CFD_FEM_QUADRANGLE_LINEAR_HPP__
#define __CFD_FEM_QUADRANGLE_LINEAR_HPP__

#include "cfd24/fem/fem_element.hpp"

namespace cfd{

///////////////////////////////////////////////////////////////////////////////
// Geometry
///////////////////////////////////////////////////////////////////////////////
class QuadrangleLinearGeometry: public IElementGeometry{
public:
	QuadrangleLinearGeometry(Point p0, Point p1, Point p2, Point p3);

	JacobiMatrix jacobi(Point xi) const override;
	Point to_physical(Point xi) const override;
private:
	Point _p0, _p1, _p2, _p3;
	Point _c_1, _c_xi, _c_eta, _c_xieta;
};

///////////////////////////////////////////////////////////////////////////////
// Basis
///////////////////////////////////////////////////////////////////////////////
class QuadrangleLinearBasis: public IElementBasis{
public:
	size_t size() const override;
	std::vector<Point> parametric_reference_points() const override;
	std::vector<BasisType> basis_types() const override;
	std::vector<double> value(Point xi) const override;
	std::vector<Vector> grad(Point xi) const override;
};

}

#endif
