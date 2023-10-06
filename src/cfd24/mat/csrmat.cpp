#include "csrmat.hpp"

using namespace cfd;

void CsrStencil::set_stencil(std::vector<size_t>&& addr, std::vector<size_t>&& cols){
	_addr = std::move(addr);
	_cols = std::move(cols);
}

size_t CsrStencil::n_nonzeros() const{
	return _cols.size();
}

size_t CsrStencil::n_rows() const{
	return _addr.size()-1;
}

const std::vector<size_t>& CsrStencil::addr() const{
	return _addr;
}

const std::vector<size_t>& CsrStencil::cols() const{
	return _cols;
}

void CsrStencil::validate() const{
	// sizes
	if (_addr.size() < 1){
		throw std::runtime_error("addr array should have more then zero entries");
	}
	if (_cols.size() != _addr.back()){
		throw std::runtime_error("cols array size should match last addr entry");
	}
	// non-decreasing
	for (size_t i=1; i<_addr.size(); ++i){
		if (_addr[i] < _addr[i-1]){
			throw std::runtime_error("addr array should be non-decreasing");
		}
	}
}

bool CsrStencil::is_in_stencil(size_t irow, size_t icol) const{
	size_t start = _addr.at(irow);
	size_t end = _addr.at(irow+1);
	for (size_t i=start; i<end; ++i){
		if (_cols.at(i) == icol){
			return true;
		}
	}
	return false;
}

double CsrStencil::value(size_t irow, size_t icol) const{
	throw std::runtime_error("CsrStencil has no values");
}

void CsrMatrix::set_values(std::vector<double>&& vals){
	_vals = std::move(vals);
}

const std::vector<double>& CsrMatrix::vals() const{
	return _vals;
}

void CsrMatrix::validate() const{
	CsrStencil::validate();

	// values size
	if (_vals.size() != n_nonzeros()){
		throw std::runtime_error("values array should have same size as the columns arrays");
	}
}

double CsrMatrix::value(size_t irow, size_t icol) const{
	const std::vector<size_t>& a = addr();
	const std::vector<size_t>& c = cols();
	size_t start = a.at(irow);
	size_t end = a.at(irow+1);
	for (size_t i=start; i<end; ++i){
		if (c.at(i) == icol){
			return _vals[i];
		}
	}
	return 0;
}

std::vector<double> CsrMatrix::mult_vec(const std::vector<double>& u) const{
	const std::vector<size_t>& a = addr();
	const std::vector<size_t>& c = cols();
	const std::vector<double>& v = vals();

	std::vector<double> ret(n_rows(), 0);
	for (size_t irow=0; irow<n_rows(); ++irow){
		size_t start = a[irow];
		size_t end = a[irow+1];
		for (size_t i=start; i<end; ++i){
			ret[irow] += v[i] * u[c[i]];
		}
	}

	return ret;
}

double CsrMatrix::mult_vec(size_t irow, const std::vector<double>& u) const{
	const std::vector<size_t>& a = addr();
	const std::vector<size_t>& c = cols();
	const std::vector<double>& v = vals();

	double ret = 0;
	size_t start = a.at(irow);
	size_t end = a.at(irow+1);
	for (size_t i=start; i<end; ++i){
		ret += v[i] * u[c[i]];
	}
	return ret;
}
