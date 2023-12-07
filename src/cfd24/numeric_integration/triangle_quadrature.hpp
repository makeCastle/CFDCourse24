#ifndef __CFD_NUMERIC_INTEGRATION_TRIANGLE_QUADRATURE_HPP__
#define __CFD_NUMERIC_INTEGRATION_TRIANGLE_QUADRATURE_HPP__

#include "cfd24/numeric_integration/quadrature.hpp"

namespace cfd{

const Quadrature* quadrature_triangle_gauss1();
const Quadrature* quadrature_triangle_gauss2();
const Quadrature* quadrature_triangle_gauss3();
const Quadrature* quadrature_triangle_gauss4();

}
#endif

