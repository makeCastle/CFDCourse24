#ifndef MAT_CSRMAT_HPP
#define MAT_CSRMAT_HPP

#include "cfd24/mat/i_sparse_mat.hpp"

namespace cfd{

/**
 * Compressed sparsed row stencil.
 * A Matrix with filled stencil and zeros values;
 */
class CsrStencil: public ISparseMatrix{
public:
	/**
	 * @brief Fills CSR stencil by address and column vectors
	 *
	 * @param addr  address vector
	 * @param cols  column vector
	 */
	void set_stencil(std::vector<size_t>&& addr, std::vector<size_t>&& cols);

	/**
	 * @brief Address array getter
	 */
	const std::vector<size_t>& addr() const;

	/**
	 * @brief Columns array getter
	 */
	const std::vector<size_t>& cols() const;

	/**
	 * Consistency check. Throws if invalid
	 */
	virtual void validate() const;

	// overrides
	size_t n_rows() const override;
	size_t n_nonzeros() const override;
	bool is_in_stencil(size_t irow, size_t icol) const override;
	double value(size_t irow, size_t icol) const override;
private:
	std::vector<size_t> _addr = {0};
	std::vector<size_t> _cols;
};

/**
 * Compressed sparsed row matrix format
 */
class CsrMatrix: public CsrStencil{
public:
	/**
	 * @brief Set matrix values
	 * 
	 * @param vals  values array
	 */
	void set_values(std::vector<double>&& vals);

	/**
	 * @brief Values array getter
	 */
	const std::vector<double>& vals() const;

	// overrides
	void validate() const override;
	double value(size_t irow, size_t icol) const override;
private:
	std::vector<double> _vals;
};


}

#endif
