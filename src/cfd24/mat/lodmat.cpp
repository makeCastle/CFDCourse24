#include "lodmat.hpp"

using namespace cfd;

LodMatrix::LodMatrix(size_t n_rows): _data(n_rows){}

const std::map<size_t, double>& LodMatrix::row(size_t i_row){
	return _data.at(i_row);
}

void LodMatrix::add_to_entry(size_t irow, size_t icol, double value){
	std::map<size_t, double>& r = _data.at(irow);
	auto found = r.find(icol);
	if (found == r.end()){
		r.emplace(icol, value);
	} else {
		found->second += value;
	}
}

void LodMatrix::remove_row(size_t irow){
	_data.at(irow).clear();
}

size_t LodMatrix::n_rows() const{
	return _data.size();
}

double LodMatrix::value(size_t irow, size_t icol) const{
	const std::map<size_t, double>& r = _data.at(irow);
	auto found = r.find(icol);
	if (found == r.end()){
		return 0;
	} else {
		return found->second;
	}
}

size_t LodMatrix::n_nonzeros() const{
	size_t ret = 0;
	for (const auto& it: _data){
		ret += it.size();
	}
	return ret;
}

bool LodMatrix::is_in_stencil(size_t irow, size_t icol) const{
	const std::map<size_t, double>& r = _data.at(irow);
	return r.find(icol) != r.end();
}
