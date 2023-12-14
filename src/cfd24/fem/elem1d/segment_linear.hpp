#ifndef __CFD_FEM_SEGMENT_LINEAR_HPP__
#define __CFD_FEM_SEGMENT_LINEAR_HPP__

#include "cfd24/fem/fem_element.hpp"

namespace cfd{

///////////////////////////////////////////////////////////////////////////////
// Geometry
///////////////////////////////////////////////////////////////////////////////
class SegmentLinearGeometry: public IElementGeometry{
public:
	SegmentLinearGeometry(Point p0, Point p1);

	JacobiMatrix jacobi(Point xi) const override;
	Point to_physical(Point xi) const override;
private:
	Point _p0, _p1;
	JacobiMatrix _jac;
};

///////////////////////////////////////////////////////////////////////////////
// Basis
///////////////////////////////////////////////////////////////////////////////
class SegmentLinearBasis: public IElementBasis{
public:
	size_t size() const override;
	std::vector<Point> parametric_reference_points() const override;
	std::vector<BasisType> basis_types() const override;
	std::vector<double> value(Point xi) const override;
	std::vector<Vector> grad(Point xi) const override;
};

///////////////////////////////////////////////////////////////////////////////
// Integrals
///////////////////////////////////////////////////////////////////////////////
class SegmentLinearIntegrals: public IElementIntegrals{
public:
	SegmentLinearIntegrals(JacobiMatrix jac);
	std::vector<double> mass_matrix() const override;
	std::vector<double> load_vector() const override;
	std::vector<double> stiff_matrix() const override;
private:
	const double _len;
};

}
#endif

