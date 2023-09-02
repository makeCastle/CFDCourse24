#include "grid1d.hpp"
#include "cfd24/grid/vtk.hpp"

using namespace cfd;

Grid1D::Grid1D(double left, double right, size_t n_cells){
	_points.push_back(Point(left));
	double h = (right - left)/n_cells;
	for (size_t i=0; i<n_cells; ++i){
		_points.push_back(Point(_points.back()[0] + h));
	}
}

size_t Grid1D::n_points() const{
	return _points.size();
}

size_t Grid1D::n_cells() const{
	return _points.size() - 1;
}

size_t Grid1D::n_faces() const{
	return _points.size();
}

Point Grid1D::point(size_t ipoint) const{
	return _points.at(ipoint);
}

std::vector<Point> Grid1D::points() const{
	return _points;
}

std::vector<size_t> Grid1D::tab_cell_point(size_t icell) const{
	return {icell, icell+1};
}

void Grid1D::save_vtk(std::string fname) const{
	std::ofstream fs(fname);
	VtkUtils::append_header("Grid1", fs);
	VtkUtils::append_points(this->points(), fs);
	//Cells
	fs << "CELLS  " << this->n_cells() << "   " << 3 * this->n_cells() << std::endl;
	for (size_t i = 0; i < this->n_cells(); ++i){
		auto v = this->tab_cell_point(i);
		fs << 2 << " " << v[0] << " " << v[1] << std::endl;
	}
	fs << "CELL_TYPES  " << this->n_cells() << std::endl;
	for (size_t i = 0; i < this->n_cells(); ++i)
		fs << 3 << std::endl;

}
