#ifndef CFD_PRINTER_HPP
#define CFD_PRINTER_HPP

#include "cfd24/mat/i_sparse_mat.hpp"

namespace cfd{

/// @brief debug procedures
namespace dbg{

/**
 * @brief prints sparse matrix to std::cout
 *
 * @param mat  sparse matrix
 *
 * prints '*' for entries that do not present in the stencil.
 */
void print(const ISparseMatrix& mat);

/**
 * @brief prints dense vector to std::cout
 *
 * @param vec  vector to print
 */
void print(const std::vector<double>& vec);

}
}
#endif
