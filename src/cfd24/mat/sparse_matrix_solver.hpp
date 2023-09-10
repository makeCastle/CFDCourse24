#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP

#include "cfd24/mat/csrmat.hpp"

namespace cfd{

/**
 * @brief Amgcl Solver for csr matrices
 */
class AmgcMatrixSolver{
public:
	/**
	 * @param maxit  maximum allowed number of iterations
	 * @param eps    tolerance value
	 */
	AmgcMatrixSolver(int maxit=1000, double eps=1e-8);
	~AmgcMatrixSolver();

	/**
	 * @brief Sets target matrix
	 *
	 * @param mat  target matrix
	 *
	 * Matrix should not be deleted until this class remains actual
	 */
	void set_matrix(const CsrMatrix& mat);

	/**
	 * @brief Sets target matrix
	 *
	 * @param mat_stencil  csr matrix stencil
	 * @param mat_values   matrix values
	 *
	 * mat_stencil and mat_values should not be deleted until this class remains actual
	 */
	void set_matrix(const CsrStencil& mat_stencil, const std::vector<double>& mat_values);

	/**
	 * @brief Solves slae Ax = rhs
	 *
	 * @param      rhs  right hand side vector
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
