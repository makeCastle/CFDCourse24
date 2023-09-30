#ifndef MAT_I_MAT_HPP
#define MAT_I_MAT_HPP

#include "cfd24/cfd_common.hpp"

namespace cfd{

/**
 * @brief Abstract matrix interface
 */
class IMatrix{
public:
	virtual ~IMatrix() = default;

	/// @brief number of rows
	virtual size_t n_rows() const = 0;

	/// @brief gets value at given address
	virtual double value(size_t irow, size_t icol) const = 0;


	/// @brief computes diagonal vector
	std::vector<double> diagonal() const{
		std::vector<double> ret(n_rows(), 0);
		for (size_t i=0; i<n_rows(); ++i){
			ret[i] = value(i, i);
		}
		return ret;
	};
};


}

#endif
