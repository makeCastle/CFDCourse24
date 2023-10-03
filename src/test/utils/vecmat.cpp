#include "vecmat.hpp"

double compute_residual(const cfd::CsrMatrix& m, const std::vector<double>& rhs, const std::vector<double>& u){
	std::vector<double> mvec = m.mult_vec(u);
	double nrm = 0;
	for (size_t irow=0; irow<m.n_rows(); ++irow){
		nrm = std::max(nrm, std::abs(mvec[irow] - rhs[irow]));
	}
	return nrm;
}

std::vector<double> vector_sum(
		const std::vector<double>& v1,
		double coef,
		const std::vector<double>& v2){

	std::vector<double> ret(v1);
	for (size_t i=0; i<ret.size(); ++i){
		ret[i] += coef * v2[i];
	}
	return ret;
}
