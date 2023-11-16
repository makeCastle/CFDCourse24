#include "searcher.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace cfd;
namespace bg = boost::geometry;

template<size_t Dim>
struct PointSearcher<Dim>::Impl{
	using point_t = bg::model::point<double, Dim, bg::cs::cartesian>;
	using value_t = std::pair<point_t, size_t>;

	void add(const std::vector<Point>& points){
		for (size_t i=0; i<points.size(); ++i){
			double x = points[i].x();
			double y = points[i].y();
			double z = points[i].z();
			_rtree.insert(std::make_pair(point_t(x, y, z), i));
		}
	}

	std::vector<size_t> nearest(const Point& p, size_t n) const{
		std::vector<value_t> vals;
		point_t point(p.x(), p.y(), p.z());
		_rtree.query(bg::index::nearest(point, n), std::back_inserter(vals));
		std::vector<size_t> ret;
		for (const auto& v: vals){
			ret.push_back(v.second);
		}
		return ret;
	}

private:
	bg::index::rtree<value_t, bg::index::quadratic<16> > _rtree;
};


template<size_t Dim>
PointSearcher<Dim>::PointSearcher(){
	_pimpl = std::make_shared<Impl>();
}

template<size_t Dim>
PointSearcher<Dim>::PointSearcher(const std::vector<Point>& points){
	_pimpl = std::make_shared<Impl>();
	add_points(points);
}

template<size_t Dim>
void PointSearcher<Dim>::add_points(const std::vector<Point>& points){
	_pimpl->add(points);
}

template<size_t Dim>
std::vector<size_t> PointSearcher<Dim>::nearest(const Point& p, size_t n) const{
	return _pimpl->nearest(p, n);
}

template struct PointSearcher<1>;
template struct PointSearcher<2>;
template struct PointSearcher<3>;
