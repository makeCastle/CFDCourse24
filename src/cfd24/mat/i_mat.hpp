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

	/// number of rows
	virtual size_t n_rows() const = 0;

	/// gets value at given address
	virtual double value(size_t irow, size_t icol) const = 0;
};


}

#endif
