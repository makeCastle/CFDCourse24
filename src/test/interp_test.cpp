#include "cfd24_test.hpp"
#include "cfd24/fem/elem1d/segment_linear.hpp"
#include "cfd24/fem/elem1d/segment_quadratic.hpp"
#include "cfd24/fem/elem1d/segment_cubic.hpp"

namespace {
struct IFunc{
	virtual double value(double) const = 0;
	virtual double dx(double) const = 0;
};

struct Func: public IFunc{
	double value(double x) const override{
		return sin(10*x*x);
	}

	double dx(double x) const override{
		return 20*x*cos(10*x*x);
	}
};

struct IInterpolator{
	IInterpolator(size_t n): _n(n), _h(1.0/n){
	}

	double value(double x) const{
		size_t icell = size_t(x/_h);
		double cell_x = 2*(x - icell*_h)/_h - 1;
		return cell_value(icell, cell_x);
	}

	double dx(double x) const{
		size_t icell = size_t(x/_h);
		double cell_x = 2*(x - icell*_h)/_h - 1;
		return cell_dx(icell, cell_x)*2.0/_h;
	}

	double point(size_t ipoint){
		return ipoint*_h;
	}

	size_t n_points() const{
		return _n + 1;
	}

	virtual double cell_value(size_t icell, double x) const = 0;
	virtual double cell_dx(size_t icell, double x) const = 0;
private:
	const size_t _n;
	const double _h;
};

struct LinearInterpolator: public IInterpolator{
	LinearInterpolator(size_t n, const IFunc& func): IInterpolator(n){
		for (size_t i=0; i<n+1; ++i){
			_vals.push_back(func.value(point(i)));
		}
	};

	double cell_value(size_t icell, double x) const override{
		std::vector<double> phi = bas.value(x);

		return _vals[icell] * phi[0] + _vals[icell+1] * phi[1];
	}

	double cell_dx(size_t icell, double x) const override{
		auto phi = bas.grad(x);

		return _vals[icell] * phi[0].x() + _vals[icell+1] * phi[1].x();
	}
private:
	std::vector<double> _vals;
	cfd::SegmentLinearBasis bas;
};

struct QuadraticInterpolator: public IInterpolator{
	QuadraticInterpolator(size_t n, const IFunc& func): IInterpolator(n){
		for (size_t i=0; i<n+1; ++i){
			_vals.push_back(func.value(point(i)));
		}
		for (size_t i=0; i<n; ++i){
			double x = (point(i) + point(i+1))/2.0;
			_vals_c.push_back(func.value(x));
		}
	};

	double cell_value(size_t icell, double x) const override{
		auto phi = bas.value(x);
		return _vals[icell] * phi[0] + _vals[icell+1] * phi[1] + _vals_c[icell] * phi[2];
	}

	double cell_dx(size_t icell, double x) const override{
		auto phi = bas.grad(x);
		return _vals[icell] * phi[0].x() + _vals[icell+1] * phi[1].x() + _vals_c[icell] * phi[2].x();
	}
private:
	std::vector<double> _vals;
	std::vector<double> _vals_c;
	cfd::SegmentQuadraticBasis bas;
};

struct CubicInterpolator: public IInterpolator{
	CubicInterpolator(size_t n, const IFunc& func): IInterpolator(n){
		for (size_t i=0; i<n+1; ++i){
			_vals.push_back(func.value(point(i)));
		}
		for (size_t i=0; i<n; ++i){
			double x0 = 2.0/3.0 * point(i) + 1.0/3.0*point(i+1);
			double x1 = 1.0/3.0 * point(i) + 2.0/3.0*point(i+1);
			_vals_c0.push_back(func.value(x0));
			_vals_c1.push_back(func.value(x1));
		}
	};

	double cell_value(size_t icell, double x) const override{
		auto phi = bas.value(x);
		return _vals[icell] * phi[0] + _vals[icell+1] * phi[1] + _vals_c0[icell] * phi[2] + _vals_c1[icell] * phi[3];
	}

	double cell_dx(size_t icell, double x) const override{
		auto phi = bas.grad(x);
		return _vals[icell] * phi[0].x() + _vals[icell+1] * phi[1].x() + _vals_c0[icell] * phi[2].x() + _vals_c1[icell] * phi[3].x();
	}
private:
	std::vector<double> _vals;
	std::vector<double> _vals_c0, _vals_c1;
	cfd::SegmentCubicBasis bas;
};

struct HermiteInterpolator: public IInterpolator{
	HermiteInterpolator(size_t n, const IFunc& func):
			IInterpolator(n),
			bas(std::make_shared<cfd::SegmentLinearGeometry>(cfd::Point{0}, cfd::Point{1.0/n})){

		for (size_t i=0; i<n+1; ++i){
			_vals.push_back(func.value(point(i)));
		}
		for (size_t i=0; i<n+1; ++i){
			_vals_dx.push_back(func.dx(point(i)));
		}
	};

	double cell_value(size_t icell, double x) const override{
		auto phi = bas.value(x);
		return _vals[icell] * phi[0] + _vals[icell+1] * phi[1] + _vals_dx[icell] * phi[2] + _vals_dx[icell+1] * phi[3];
	}

	double cell_dx(size_t icell, double x) const override{
		auto phi = bas.grad(x);
		return _vals[icell] * phi[0].x() + _vals[icell+1] * phi[1].x() + _vals_dx[icell] * phi[2].x() + _vals_dx[icell+1] * phi[3].x();
	}
private:
	std::vector<double> _vals;
	std::vector<double> _vals_dx;
	cfd::SegmentHermiteBasis bas;
};

template<typename TInterp>
void print_norm(const IFunc& func){
	constexpr size_t ntest = 100'000;
	constexpr double htest = 1.0/ntest;

	for (size_t n: {10, 20, 50, 100, 200, 500, 1000}){
		TInterp interp(n, func);
		double diff = 0;
		for (size_t i=0; i<ntest; ++i){
			double x = i*htest;
			double d = func.value(x) - interp.value(x);
			double coef = (i == 0 || i == ntest-1) ? 0.5 : 1.0;
			diff += coef * htest * d*d;
		}
		diff = std::sqrt(diff);
		std::cout << n << " " << diff << std::endl;
	}
};

template<typename TInterp>
void print_dx_norm(const IFunc& func){
	constexpr size_t ntest = 10'000;
	constexpr double htest = 1.0/ntest;

	for (size_t n: {10, 20, 50, 100, 200, 500, 1000}){
		TInterp interp(n, func);
		double diff = 0;
		for (size_t i=0; i<ntest; ++i){
			double x = i*htest;
			double d = func.dx(x) - interp.dx(x);
			double coef = (i == 0 || i == ntest-1) ? 0.5 : 1.0;
			diff += coef * htest * d * d;
		}
		diff = std::sqrt(diff);
		std::cout << n << " " << diff << std::endl;
	}
};
}

TEST_CASE("MY", "[.my]"){
	Func func;

	// =============== func norm
	std::cout << "LINEAR" << std::endl;
	print_norm<LinearInterpolator>(func);

	std::cout << "QUADRATIC" << std::endl;
	print_norm<QuadraticInterpolator>(func);

	std::cout << "CUBIC" << std::endl;
	print_norm<CubicInterpolator>(func);

	std::cout << "HERMITE" << std::endl;
	print_norm<HermiteInterpolator>(func);

	// =============== dx norm
	std::cout << "LINEAR" << std::endl;
	print_dx_norm<LinearInterpolator>(func);

	std::cout << "QUADRATIC" << std::endl;
	print_dx_norm<QuadraticInterpolator>(func);

	std::cout << "CUBIC" << std::endl;
	print_dx_norm<CubicInterpolator>(func);

	std::cout << "HERMITE" << std::endl;
	print_dx_norm<HermiteInterpolator>(func);
}
