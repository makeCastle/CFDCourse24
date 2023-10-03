#include "regular_grid2d.hpp"
#include "cfd24/grid/vtk.hpp"

using namespace cfd;

RegularGrid2D::RegularGrid2D(double x0, double x1, double y0, double y1, size_t nx, size_t ny){
	_x.push_back(x0);
	double hx = (x1 - x0)/nx;
	for (size_t i=0; i<nx; ++i){
		_x.push_back(_x.back() + hx);
	}
	_y.push_back(y0);
	double hy = (y1 - y0)/ny;
	for (size_t i=0; i<ny; ++i){
		_y.push_back(_y.back() + hy);
	}
}

RegularGrid2D::RegularGrid2D(const std::vector<double>& x, const std::vector<double>& y):
	_x(x), _y(y){ }

size_t RegularGrid2D::nx() const{
	return _x.size() - 1;
}

size_t RegularGrid2D::ny() const{
	return _y.size() - 1;
}

size_t RegularGrid2D::n_points() const{
	return _x.size() * _y.size();
}

size_t RegularGrid2D::n_cells() const{
	return (_x.size()-1)*(_y.size()-1);
}

size_t RegularGrid2D::n_faces() const{
	_THROW_NOT_IMP_;
}

Point RegularGrid2D::point(size_t ipoint) const{
	size_t ix = ipoint % _x.size();
	size_t iy = ipoint / _x.size();
	return Point(_x[ix], _y[iy]);
}

std::vector<Point> RegularGrid2D::points() const{
	std::vector<Point> ret;
	for (size_t i=0; i<n_points(); ++i){
		ret.push_back(point(i));
	}
	return ret;
}

std::vector<size_t> RegularGrid2D::tab_cell_point(size_t icell) const{
	size_t ix = icell % (_x.size() - 1);
	size_t iy = icell / (_x.size() - 1);

	size_t n0 = ix + iy * _x.size();
	size_t n1 = ix + 1 + iy * _x.size();
	size_t n2 = ix + 1 + (iy + 1) * _x.size();
	size_t n3 = ix + (iy + 1) * _x.size();
	return {n0, n1, n2, n3};
}

void RegularGrid2D::save_vtk(std::string fname) const{
	std::ofstream fs(fname);
	VtkUtils::append_header("Grid2", fs);
	VtkUtils::append_points(this->points(), fs);
	//Cells
	fs << "CELLS  " << this->n_cells() << "   " << 5 * this->n_cells() << std::endl;
	for (size_t i = 0; i < this->n_cells(); ++i){
		auto v = this->tab_cell_point(i);
		fs << 4 << " " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;
	}
	fs << "CELL_TYPES  " << this->n_cells() << std::endl;
	for (size_t i = 0; i < this->n_cells(); ++i)
		fs << 9 << std::endl;

}

size_t RegularGrid2D::to_linear_point_index(split_index_t point_split_index) const{
	return point_split_index[0] + _x.size() * point_split_index[1];
}

auto RegularGrid2D::to_split_point_index(size_t ipoint) const -> split_index_t{
	return split_index_t{
		ipoint % _x.size(),
		ipoint / _x.size()
	};
}

size_t RegularGrid2D::cell_centered_grid_index_ip_jp(size_t i, size_t j) const{
	if (i >= _x.size()-1 || j >= _y.size()-1){
		std::ostringstream oss;
		oss << "Invalid RegularGrid2d::cell_centered_grid_index_ip_jp() arguments: ";
		oss << std::endl << "    ";
		oss << "i=" << i << ", j=" << j;
		throw std::runtime_error(oss.str());
	}
	return i + j*(_x.size() - 1);
}

size_t RegularGrid2D::xface_grid_index_ip_j(size_t i, size_t j) const{
	if (i >= _x.size()-1 || j >= _y.size()){
		std::ostringstream oss;
		oss << "Invalid RegularGrid2d::xface_grid_index_ip_j() arguments: ";
		oss << std::endl << "    ";
		oss << "i=" << i << ", j=" << j;
		throw std::runtime_error(oss.str());
	}
	return i + j*(_x.size() - 1);
}

size_t RegularGrid2D::yface_grid_index_i_jp(size_t i, size_t j) const{
	if (i >= _x.size() || j >= _y.size() - 1){
		std::ostringstream oss;
		oss << "Invalid RegularGrid2d::yface_grid_index_ip_j() arguments: ";
		oss << std::endl << "    ";
		oss << "i=" << i << ", j=" << j;
		throw std::runtime_error(oss.str());
	}
	return i + j*_x.size();
}

RegularGrid2D RegularGrid2D::cell_centered_grid() const{
	std::vector<double> ret_x;
	for (size_t i=0; i<_x.size()-1; ++i){
		ret_x.push_back((_x[i] + _x[i+1])/2);
	}
	std::vector<double> ret_y;
	for (size_t i=0; i<_y.size()-1; ++i){
		ret_y.push_back((_y[i] + _y[i+1])/2);
	}
	return RegularGrid2D(ret_x, ret_y);
}

RegularGrid2D RegularGrid2D::xface_centered_grid() const{
	std::vector<double> ret_x;
	for (size_t i=0; i<_x.size()-1; ++i){
		ret_x.push_back((_x[i] + _x[i+1])/2);
	}
	return RegularGrid2D(ret_x, _y);
}

RegularGrid2D RegularGrid2D::yface_centered_grid() const{
	std::vector<double> ret_y;
	for (size_t i=0; i<_y.size()-1; ++i){
		ret_y.push_back((_y[i] + _y[i+1])/2);
	}
	return RegularGrid2D(_x, ret_y);
}
