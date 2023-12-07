#ifndef __CFD_FEM_I_ELEMENT_HPP__
#define __CFD_FEM_I_ELEMENT_HPP__

#include "cfd24/mat/densemat.hpp"
#include "cfd24/grid/i_grid.hpp"
#include "cfd24/geom/jacobi.hpp"

namespace cfd{

///////////////////////////////////////////////////////////////////////////////
// Element Geometry
///////////////////////////////////////////////////////////////////////////////
class IElementGeometry{
public:
	virtual ~IElementGeometry() = default;

	virtual JacobiMatrix jacobi(Point xi) const = 0;
	virtual Point to_physical(Point xi) const  { _THROW_NOT_IMP_; }
	virtual Point to_parametric(Point p) const { _THROW_NOT_IMP_; }
};

///////////////////////////////////////////////////////////////////////////////
// Element Basis
///////////////////////////////////////////////////////////////////////////////
class IElementBasis{
public:
	virtual ~IElementBasis() = default;

	virtual size_t size() const = 0;
	virtual std::vector<Point> parametric_reference_points() const = 0;
	virtual std::vector<double> value(Point xi) const = 0;
	virtual std::vector<Vector> grad(Point xi) const = 0;
};

///////////////////////////////////////////////////////////////////////////////
// Element Integrals
///////////////////////////////////////////////////////////////////////////////
class IElementIntegrals{
public:
	virtual ~IElementIntegrals() = default;

	virtual std::vector<double> load_vector() const  { _THROW_NOT_IMP_; }
	virtual std::vector<double> mass_matrix() const  { _THROW_NOT_IMP_; }
	virtual std::vector<double> stiff_matrix() const { _THROW_NOT_IMP_; }
};

///////////////////////////////////////////////////////////////////////////////
// FemElement
///////////////////////////////////////////////////////////////////////////////
struct FemElement{
	std::shared_ptr<const IElementGeometry> geometry;
	std::shared_ptr<const IElementBasis> basis;
	std::shared_ptr<const IElementIntegrals> integrals;
};

}
#endif
