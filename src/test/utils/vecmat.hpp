#ifndef __VECMAT_HPP__
#define __VECMAT_HPP__

#include <vector>
#include "cfd24/mat/csrmat.hpp"

// -> max(abs(m*u - rhs))
double compute_residual(const cfd::CsrMatrix& m, const std::vector<double>& rhs, const std::vector<double>& u);

// -> v1 + coef*v2
std::vector<double> vector_sum(
		const std::vector<double>& v1,
		double coef,
		const std::vector<double>& v2);

#endif
