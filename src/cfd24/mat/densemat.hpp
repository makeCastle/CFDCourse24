#ifndef CFD_DENSE_MAT_HPP
#define CFD_DENSE_MAT_HPP

#include "cfd24/mat/i_mat.hpp"

namespace cfd{

class DenseMatrix: public IMatrix{
public:
	DenseMatrix(size_t nrows, size_t ncols);
	DenseMatrix(size_t nrows, size_t ncols, const std::vector<double>& values);
	
	void set_value(size_t irow, size_t icol, double value);

	DenseMatrix transpose() const;
	DenseMatrix mult_mat(const DenseMatrix& mat) const;
	DenseMatrix inverse() const;

	size_t n_cols() const;

	// overriden
	size_t n_rows() const override;
	double value(size_t irow, size_t icol) const override;
	std::vector<double> mult_vec(const std::vector<double>& u) const override;
	double mult_vec(size_t irow, const std::vector<double>& u) const override;
private:
	const size_t _nrows;
	const size_t _ncols;
	std::vector<double> _data;
};

}
#endif
