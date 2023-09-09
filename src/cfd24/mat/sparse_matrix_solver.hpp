#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP

#include "cfd24/mat/csrmat.hpp"

namespace cfd{

/**
 * Amgcl Solver for csr matrices
 */
class AmgcMatrixSolver{
public:
	/**
	 * Constructor
	 *
	 * @param maxit  maximum number of iterations
	 * @param eps    tolerance value
	 */
	AmgcMatrixSolver(int maxit=1000, double eps=1e-8);
	~AmgcMatrixSolver();

	/**
	 * Sets target matrix
	 *
	 * @param mat  target matrix
	 *
	 * Matrix should not be deleted until this class remains actual
	 */
	void set_matrix(const CsrMatrix& mat);

	/**
	 * Sets target matrix
	 *
	 * @param mat_stencil  csr matrix stencil
	 * @param mat_value    matrix values
	 *
	 * mat_stencil and mat_value should not be deleted until this class remains actual
	 */
	void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values);

	/**
	 * Solves slae Ax = rhs
	 *
	 * @param rhs  right hand side vector
	 * @param[out] ret  output vector
	 *
	 * Given values of output vector will be used as the initial values.
	 * If it is empty, solution will be initialized with zeros.
	 */
	void solve(const std::vector<double>& rhs, std::vector<double>& ret) const;
private:
	int _maxit;
	double _tolerance;

	struct Impl;
	std::unique_ptr<Impl> _pimpl;
};

}

#endif
