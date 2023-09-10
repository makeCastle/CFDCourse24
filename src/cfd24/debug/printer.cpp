#include "printer.hpp"
#include <iomanip>

using namespace cfd;

void dbg::print(const ISparseMatrix& mat){
	size_t n = mat.n_rows();
	std::cout << "-- SIZE = " << n << "x" << n <<std::endl;
	for (size_t irow=0; irow < n; ++irow){
		for (size_t icol=0; icol < n; ++icol){
			std::cout << std::setw(6);
			if (mat.is_in_stencil(irow, icol)){
				std::cout << mat.value(irow, icol);
			} else {
				std::cout << "*";
			}
			std::cout << " ";
		}
		std::cout << std::endl;
	}
}

void dbg::print(const std::vector<double>& vec){
	std::cout << "-- SIZE = " << vec.size() << std::endl;
	for (double v: vec){
		std::cout << v << std::endl;
	}
}
