#ifndef MAT_LODMAT_HPP
#define MAT_LODMAT_HPP

#include "cfd24/mat/i_sparse_mat.hpp"

namespace cfd{

/**
 * List of dictionaries sparse matrix format
 */
class LodMatrix: public ISparseMatrix{
public:
	/**
	 * Construct empty matrix
	 *
	 * @param n_rows number of rows
	 */
	LodMatrix(size_t n_rows);

	/**
	 * Get row values
	 *
	 * @param irow  row index
	 * @return column->value dictionary
	 */
	const std::map<size_t, double>& row(size_t irow);

	/**
	 * Adds value to the given matrix entry
	 *
	 * @param irow   row index
	 * @param icol   column index
	 * @param value  value to be added
	 *
	 * Performs Matrix[i, j] += value operation
	 */
	void add_to_entry(size_t irow, size_t icol, double value);

	/**
	 * Removes row from sparse matrix stencil
	 *
	 * @param irow  row index
	 *
	 * All values at the given row will be equal to zero
	 */
	void remove_row(size_t irow);

	// overrides
	size_t n_rows() const override;
	double value(size_t irow, size_t icol) const override;
	size_t n_nonzeros() const override;
	bool is_in_stencil(size_t irow, size_t icol) const override;
private:
	std::vector<std::map<size_t, double>> _data;
};

}

#endif
