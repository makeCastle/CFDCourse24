#include "cfd24_test.hpp"
#include "cfd24/grid/grid1d.hpp"
#include "cfd24/mat/csrmat.hpp"
#include "cfd24/mat/sparse_matrix_solver.hpp"
#include "cfd24/grid/vtk.hpp"
#include "cfd24/mat/lodmat.hpp"
#include <cfd24/grid/regular_grid2d.hpp>


using namespace cfd;

class ATestTransport1Worker {
public:
	static double init_solution(double x) {
		//return (x >= -0.1 && x <= 0.1) ? 1.0 : 0.0;
		constexpr double sigma = 0.1;
		return exp(-x * x / sigma / sigma);
	}
	double exact_solution(double x) {
		return init_solution(x - _time);
	}

	ATestTransport1Worker(size_t n_cells) : _grid(0, 1, n_cells), _h(1.0 / n_cells), _u(_grid.n_points()) {
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			_u[i] = init_solution(_grid.point(i).x());
		}
	}

	double step(double tau) {
		_time += tau;
		impl_step(tau);
		return compute_norm2();
	}

	void save_vtk(const std::string& filename) {
		// save grid
		_grid.save_vtk(filename);

		// save numerical solution
		VtkUtils::add_point_data(_u, "numerical", filename);

		// save exact solution
		std::vector<double> exact(_grid.n_points());
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			exact[i] = exact_solution(_grid.point(i).x());
		}
		VtkUtils::add_point_data(exact, "exact", filename);
	}

	double current_time() const {
		return _time;
	}

	double compute_norm2() {
		// weights
		std::vector<double> w(_grid.n_points(), _h);
		w[0] = w[_grid.n_points() - 1] = _h / 2;

		// sum
		double sum = 0;
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			double diff = _u[i] - exact_solution(_grid.point(i).x());
			sum += w[i] * diff * diff;
		}

		double len = _grid.point(_grid.n_points() - 1).x() - _grid.point(0).x();
		return std::sqrt(sum / len);
	}
protected:
	Grid1D _grid;
	double _h;
	std::vector<double> _u;
	double _time = 0;

private:
	virtual void impl_step(double tau) = 0;
};

///////////////////////////////////////////////////////////////////////////////
// Explicit transport 1D solver
///////////////////////////////////////////////////////////////////////////////

class TestTransport1WorkerExplicit : public ATestTransport1Worker {
public:
	TestTransport1WorkerExplicit(size_t n_cells) : ATestTransport1Worker(n_cells) {}

private:
	void impl_step(double tau) override {
		std::vector<double> u_old(_u);
		_u[0] = exact_solution(_grid.point(0).x());
		for (size_t i = 1; i < _grid.n_points(); ++i) {
			_u[i] = u_old[i] - tau / _h * (u_old[i] - u_old[i - 1]);
		}
	}
};

TEST_CASE("Transport 1D solver, explicit", "[transport1-explicit]") {
	std::cout << std::endl << "--- cfd24_test [transport1-explicit] --- " << std::endl;
	// parameters
	const double tend = 0.5;
	const double V = 1.0;
	const double L = 1.0;
	size_t n_cells = 100;
	double Cu = 0.9;
	double h = L / n_cells;
	double tau = Cu * h / V;

	// solver
	TestTransport1WorkerExplicit worker(n_cells);

	// saver
	VtkUtils::TimeDependentWriter writer("transport1-explicit");
	std::string out_filename = writer.add(0);
	worker.save_vtk(out_filename);

	double norm = 0;
	while (worker.current_time() < tend - 1e-6) {
		// solve problem
		norm = worker.step(tau);
		// export solution to vtk
		out_filename = writer.add(worker.current_time());
		worker.save_vtk(out_filename);
	};
	std::cout << 1.0 / tau << " " << norm << std::endl;
	CHECK(norm == Approx(0.0138123932).margin(1e-5));
}

///////////////////////////////////////////////////////////////////////////////
// Implicit transport 1D solver
///////////////////////////////////////////////////////////////////////////////

class TestTransport1WorkerImplicit : public ATestTransport1Worker {
public:
	TestTransport1WorkerImplicit(size_t n_cells) : ATestTransport1Worker(n_cells) {}

private:
	void impl_step(double tau) override {
		AmgcMatrixSolver& slv = build_solver(tau);
		std::vector<double> rhs = build_rhs(tau);
		slv.solve(rhs, _u);
	}

	AmgcMatrixSolver _solver;
	double _last_used_tau = 0;

	AmgcMatrixSolver& build_solver(double tau) {
		if (_last_used_tau != tau) {
			CsrMatrix mat = build_lhs(tau);
			_solver.set_matrix(mat);
			_last_used_tau = tau;
		}
		return _solver;
	}

	virtual CsrMatrix build_lhs(double tau) {
		LodMatrix mat(_u.size());
		mat.set_value(0, 0, 1.0);
		mat.set_value(_u.size() - 1, _u.size() - 1, 1.0);
		double diag = 1.0 + tau / _h;
		double nondiag = -tau / _h;
		for (size_t i = 1; i < _u.size() - 1; ++i) {
			mat.set_value(i, i, diag);
			mat.set_value(i, i - 1, nondiag);
		}
		return mat.to_csr();
	}

	virtual std::vector<double> build_rhs(double tau) {
		std::vector<double> rhs(_u);
		rhs[0] = exact_solution(_grid.point(0).x());
		rhs.back() = exact_solution(_grid.point(_grid.n_points() - 1).x());
		return rhs;
	}
};

TEST_CASE("Transport 1D solver, implicit", "[transport1-implicit]") {
	std::cout << std::endl << "--- cfd24_test [transport1-implicit] --- " << std::endl;
	const double tend = 0.5;
	const double V = 1.0;
	const double L = 1.0;
	size_t n_cells = 100;
	double Cu = 0.9;

	double h = L / n_cells;
	double tau = Cu * h / V;
	TestTransport1WorkerImplicit worker(n_cells);
	VtkUtils::TimeDependentWriter writer("transport1-implicit");
	double norm;
	worker.save_vtk(writer.add(0));
	while (worker.current_time() < tend - 1e-6) {
		norm = worker.step(tau);
		worker.save_vtk(writer.add(worker.current_time()));
	};
	std::cout << 1.0 / tau << " " << norm << std::endl;
	CHECK(norm == Approx(0.137664).margin(1e-5));
}

///////////////////////////////////////////////////////////////////////////////
// Crank-Nicolson transport 1D solver
///////////////////////////////////////////////////////////////////////////////

class TestTransport1WorkerCN : public TestTransport1WorkerImplicit {
public:
	TestTransport1WorkerCN(size_t n_cells) : TestTransport1WorkerImplicit(n_cells) {}
private:
	virtual CsrMatrix build_lhs(double tau) {
		LodMatrix mat(_u.size());
		mat.set_value(0, 0, 1.0);
		mat.set_value(_u.size() - 1, _u.size() - 1, 1.0);
		double diag = 1.0 + 0.5 * tau / _h;
		double nondiag = -0.5 * tau / _h;
		for (size_t i = 1; i < _u.size() - 1; ++i) {
			mat.set_value(i, i, diag);
			mat.set_value(i, i - 1, nondiag);
		}
		return mat.to_csr();
	}

	virtual std::vector<double> build_rhs(double tau) {
		std::vector<double> rhs(_u);
		rhs[0] = exact_solution(_grid.point(0).x());
		rhs.back() = exact_solution(_grid.point(_grid.n_points() - 1).x());
		for (size_t i = 1; i < rhs.size() - 1; ++i) {
			rhs[i] -= 0.5 * tau / _h * (_u[i] - _u[i - 1]);
		}
		return rhs;
	}
};

TEST_CASE("Transport 1D solver, Crank-Nicolson", "[transport1-cn]") {
	std::cout << std::endl << "--- cfd24_test [transport1-cn] --- " << std::endl;
	const double tend = 0.5;
	const double V = 1.0;
	const double L = 1.0;
	size_t n_cells = 100;
	double Cu = 0.9;

	double h = L / n_cells;
	double tau = Cu * h / V;
	TestTransport1WorkerCN worker(n_cells);
	VtkUtils::TimeDependentWriter writer("transport1-cn");
	double norm;
	worker.save_vtk(writer.add(0));
	while (worker.current_time() < tend - 1e-6) {
		norm = worker.step(tau);
		worker.save_vtk(writer.add(worker.current_time()));
	};
	std::cout << 1.0 / tau << " " << norm << std::endl;
	CHECK(norm == Approx(0.0937748).margin(1e-5));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ATestTransport2Worker {
public:
	static double init_solution(double x, double y) {
		//return (x >= -0.1 && x <= 0.1) ? 1.0 : 0.0;

		double u = 0.0;
		if (x >= -1.0 && x <= -0.8 && y >= -0.1 && y <= 0.1) {
			u = 1.0;
		}

		return u;
	}
	double exact_solution(double x, double y) {
		return init_solution(x - _time, y);
	}

	ATestTransport2Worker(size_t n_cells) : _grid(-1.0, 1.0, -1.0, 1.0, n_cells, n_cells), _hx(2.0 / n_cells), _hy(2.0 / n_cells), _u(_grid.n_points()) {
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			_u[i] = init_solution(_grid.point(i).x(), _grid.point(i).y());
		}
	}

	double step(double tau, double U, double V) {
		_time += tau;
		impl_step(tau, U, V);
		return compute_norm2();
	}

	void save_vtk(const std::string& filename) {
		// save grid
		_grid.save_vtk(filename);

		// save numerical solution
		VtkUtils::add_point_data(_u, "numerical", filename);

		// save exact solution
		std::vector<double> exact(_grid.n_points());
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			exact[i] = exact_solution(_grid.point(i).x(), _grid.point(i).y());
		}
		VtkUtils::add_point_data(exact, "exact", filename);
	}

	double current_time() const {
		return _time;
	}

	double compute_norm2() {
		// weights
		std::vector<double> w(_grid.n_points(), _hx * _hy);

		unsigned long long ny = _grid.point(_grid.n_points() - 1).y() / _hy;
		unsigned long long nx = _grid.point(_grid.n_points() - 1).x() / _hx;
		w[_grid.to_linear_point_index({ 0, 0 })] = w[_grid.to_linear_point_index({ 0,  ny })] = w[_grid.to_linear_point_index({ nx,  ny })] =
			w[_grid.to_linear_point_index({ nx,  0 })] = (_hx * _hy) / 4.0;
		for (size_t i = 1; i < nx; ++i)
		{
			w[_grid.to_linear_point_index({ i, 0 })] = (_hx * _hy) / 2.0;
		}
		for (size_t i = 1; i < ny; ++i)
		{
			w[_grid.to_linear_point_index({ nx, i })] = (_hx * _hy) / 2.0;
		}
		for (size_t i = 1; i < nx; ++i)
		{
			w[_grid.to_linear_point_index({ i, ny })] = (_hx * _hy) / 2.0;
		}
		for (size_t i = 1; i < ny; ++i)
		{
			w[_grid.to_linear_point_index({ 0, i })] = (_hx * _hy) / 2.0;
		}

		// sum
		double sum = 0;
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			double diff = _u[i] - exact_solution(_grid.point(i).x(), _grid.point(i).y());
			sum += w[i] * diff * diff;
		}

		double lenx = _grid.point(_grid.n_points() - 1).x() - _grid.point(0).x();
		double leny = _grid.point(_grid.n_points() - 1).y() - _grid.point(0).y();
		return std::sqrt(sum / (lenx * leny));
	}
protected:
	RegularGrid2D _grid;
	double _hx;
	double _hy;
	std::vector<double> _u;
	double _time = 0;

private:
	virtual void impl_step(double tau, double U, double V) = 0;
};

///////////////////////////////////////////////////////////////////////////////
// Explicit transport 2D solver 1 task
///////////////////////////////////////////////////////////////////////////////
class TestTransport2WorkerExplicit : public ATestTransport2Worker {
public:
	TestTransport2WorkerExplicit(size_t n_cells) : ATestTransport2Worker(n_cells), pointsx(n_cells + 1), pointsy(n_cells + 1) {}

private:
	size_t pointsx;
	size_t pointsy;
	void impl_step(double tau, double U, double V) override {
		size_t upx = 0;
		size_t upy = 0;
		std::vector<double> u_old(_u);
		for (size_t i = 0; i < pointsy; ++i) {
			size_t ind = _grid.to_linear_point_index({ 0, i });
			_u[ind] = exact_solution(_grid.point(ind).x(), _grid.point(ind).y());
		}

		for (size_t i = 1; i < pointsx; ++i) {
			for (size_t j = 1; j < pointsy; ++j) {
				size_t k = _grid.to_linear_point_index({ i, j });
				if (U >= 0) {
					upx = _grid.to_linear_point_index({ i - 1, j });
				}
				else {
					upx = _grid.to_linear_point_index({ i + 1, j });
				}
				if (V >= 0) {
					upy = _grid.to_linear_point_index({ i, j - 1 });
				}
				else {
					upy = _grid.to_linear_point_index({ i, j + 1 });
				}
				_u[k] = u_old[k] - tau * (abs(U) * (_u[k] - _u[upx]) / _hx + abs(V) * (_u[k] - _u[upy]) / _hy);
			}

		}
	}
};

TEST_CASE("Transport 2D solver, explicit", "[transport2-explicit]") {
	std::cout << std::endl << "--- cfd24_test [transport2-explicit] --- " << std::endl;
	// parameters
	const double tend = 0.5;
	const double V = 0.0;
	const double U = 1.0;
	const double L1 = 1.0;
	const double L2 = 1.0;
	size_t n_cells = 40;
	double Cu = 0.9;
	double hx = L1 / n_cells;
	double hy = L2 / n_cells;
	//double tau = Cu * h / V; //find 

	//double tau = hx * 0.9;
	double tau = 0.05494;

	// solver
	TestTransport2WorkerExplicit worker(n_cells);

	// saver
	VtkUtils::TimeDependentWriter writer("transport2-explicit");
	std::string out_filename = writer.add(0);
	worker.save_vtk(out_filename);

	double norm = 0;
	while (worker.current_time() < tend - 1e-6) {
		// solve problem
		norm = worker.step(tau, U, V);
		// export solution to vtk
		out_filename = writer.add(worker.current_time());
		worker.save_vtk(out_filename);
		//std::cout << 1.0 / tau << " " << norm << std::endl;
	};
	std::cout << 1.0 / tau << " " << norm << std::endl;
	//CHECK(norm == Approx(0.0138123932).margin(1e-5));
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ATestTransport3Worker {
public:
	static double init_solution(double x, double y) {
		//return (x >= -0.1 && x <= 0.1) ? 1.0 : 0.0;
		//double e = 2.7182818;
		double r = 0.0;
		double u = 0.0;
		double sigma = 0.1;
		//if (x >= -1.0 && x <= -0.8 && y >= -0.1 && y <= 0.1) {
		r = sqrt((x - 0.5) * (x - 0.5) + y * y);
		u = exp(-r * r / (sigma * sigma));
		//}
		return u;
	}
	double exact_solution(double x, double y) {
		double xc = 0.5 * cos(0.5 * _time);
		double yc = 0.5 * sin(0.5 * _time);
		double r = 0.0;
		double u = 0.0;
		double sigma = 0.1;
		//	double e = 2.7182818;
			//if (x >= -1.0 && x <= -0.8 && y >= -0.1 && y <= 0.1) {
		r = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
		u = exp(-r * r / (sigma * sigma));
		//}
		return u;
	}


	ATestTransport3Worker(size_t n_cells) : _grid(-1.0, 1.0, -1.0, 1.0, n_cells, n_cells), _hx(2.0 / n_cells), _hy(2.0 / n_cells), _u(_grid.n_points()) {
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			_u[i] = init_solution(_grid.point(i).x(), _grid.point(i).y());
		}
	}

	double step(double tau) {
		_time += tau;
		impl_step(tau);
		return compute_norm2();
	}

	void save_vtk(const std::string& filename) {
		// save grid
		_grid.save_vtk(filename);

		// save numerical solution
		VtkUtils::add_point_data(_u, "numerical", filename);

		// save exact solution
		std::vector<double> exact(_grid.n_points());
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			exact[i] = exact_solution(_grid.point(i).x(), _grid.point(i).y());
		}
		VtkUtils::add_point_data(exact, "exact", filename);
	}

	double current_time() const {
		return _time;
	}

	double compute_norm2() {
		// weights
		std::vector<double> w(_grid.n_points(), _hx * _hy);

		unsigned long long ny = _grid.point(_grid.n_points() - 1).y() / _hy;
		unsigned long long nx = _grid.point(_grid.n_points() - 1).x() / _hx;
		w[_grid.to_linear_point_index({ 0, 0 })] = w[_grid.to_linear_point_index({ 0,  ny })] = w[_grid.to_linear_point_index({ nx,  ny })] =
			w[_grid.to_linear_point_index({ nx,  0 })] = (_hx * _hy) / 4.0;
		for (size_t i = 1; i < nx; ++i)
		{
			w[_grid.to_linear_point_index({ i, 0 })] = (_hx * _hy) / 2.0;
		}
		for (size_t i = 1; i < ny; ++i)
		{
			w[_grid.to_linear_point_index({ nx, i })] = (_hx * _hy) / 2.0;
		}
		for (size_t i = 1; i < nx; ++i)
		{
			w[_grid.to_linear_point_index({ i, ny })] = (_hx * _hy) / 2.0;
		}
		for (size_t i = 1; i < ny; ++i)
		{
			w[_grid.to_linear_point_index({ 0, i })] = (_hx * _hy) / 2.0;
		}

		// sum
		double sum = 0;
		for (size_t i = 0; i < _grid.n_points(); ++i) {
			double diff = _u[i] - exact_solution(_grid.point(i).x(), _grid.point(i).y());
			sum += w[i] * diff * diff;
		}

		double lenx = _grid.point(_grid.n_points() - 1).x() - _grid.point(0).x();
		double leny = _grid.point(_grid.n_points() - 1).y() - _grid.point(0).y();
		return std::sqrt(sum / (lenx * leny));
	}
protected:
	RegularGrid2D _grid;
	double _hx;
	double _hy;
	std::vector<double> _u;
	double _time = 0;

private:
	virtual void impl_step(double tau) = 0;
};
///////////////////////////////////////////////////////////////////////////////
// Explicit transport 2D solver 2task
///////////////////////////////////////////////////////////////////////////////
class TestTransport3WorkerExplicit : public ATestTransport3Worker {
public:
	TestTransport3WorkerExplicit(size_t n_cells) : ATestTransport3Worker(n_cells), pointsx(n_cells + 1), pointsy(n_cells + 1) {}

	double U0(double y) {
		double u;

		u = -y;

		return u;
	}

	double V0(double x) {
		double v;

		v = x;

		return v;
	}
	/*
	std::vector<double> U0() {
		std::vector<double> u;
		for (size_t i = 0; i < pointsx; ++i) {
			for (size_t j = 0; j < pointsy; ++j) {
				u.push_back(-_grid.point(_grid.to_linear_point_index({ i, j })).y());
			}
		}
		return u;
	}

	std::vector<double> V0() {
		std::vector<double> v;
		for (size_t i = 0; i < pointsx; ++i) {
			for (size_t j = 0; j < pointsy; ++j) {
				v.push_back(_grid.point(_grid.to_linear_point_index({ i, j })).x());
			}
		}
		return v;
	}
	*/

private:
	size_t pointsx;
	size_t pointsy;
	void impl_step(double tau) override {
		size_t upx = 0;
		size_t upy = 0;
		std::vector<double> u_old(_u);
		for (size_t i = 0; i < pointsy; ++i) {
			size_t ind = _grid.to_linear_point_index({ 0, i }); //left boundary
			_u[ind] = exact_solution(_grid.point(ind).x(), _grid.point(ind).y());
		}
		for (size_t i = 0; i < pointsx; ++i) {
			size_t ind = _grid.to_linear_point_index({ i, 0 }); //bottom boundary
			_u[ind] = exact_solution(_grid.point(ind).x(), _grid.point(ind).y());
		}
		for (size_t i = 0; i < pointsx; ++i) {
			size_t ind = _grid.to_linear_point_index({ i, pointsy - 1 }); //top boundary
			_u[ind] = exact_solution(_grid.point(ind).x(), _grid.point(ind).y());
		}
		for (size_t i = 0; i < pointsy; ++i) {
			size_t ind = _grid.to_linear_point_index({ pointsx - 1, i }); //right boundary
			_u[ind] = exact_solution(_grid.point(ind).x(), _grid.point(ind).y());
		}

		//int a = pointsx;
		//int b = pointsy;

		for (size_t i = 1; i < pointsx - 1; ++i) {
			for (size_t j = 1; j < pointsy - 1; ++j) {

				size_t k = _grid.to_linear_point_index({ i, j });
				double x = _grid.point(k).x();
				double y = _grid.point(k).y();
				double U = U0(y);
				double V = V0(x);
				if (U >= 0) {
					upx = _grid.to_linear_point_index({ i - 1, j });
				}
				else {
					upx = _grid.to_linear_point_index({ i + 1, j });
				}
				if (V >= 0) {
					upy = _grid.to_linear_point_index({ i, j - 1 });
				}
				else {
					upy = _grid.to_linear_point_index({ i, j + 1 });
				}
				_u[k] = u_old[k] - tau * (abs(U) * (_u[k] - _u[upx]) / _hx + abs(V) * (_u[k] - _u[upy]) / _hy);
			}

		}
	}
};

TEST_CASE("Transport 2D solver 2 task, explicit", "[transport3-explicit]") {
	std::cout << std::endl << "--- cfd24_test [transport3-explicit] --- " << std::endl;
	// parameters
	const double tend = 0.5; //конечное время
	std::vector<double> V;
	std::vector<double> U;
	const double L1 = 1.0;
	const double L2 = 1.0;
	size_t n_cells = 40;
	double Cu = 0.9;
	double hx = L1 / n_cells;
	double hy = L2 / n_cells;
	//double tau = Cu * h / V; //find 

	/*double tau = hx * 0.9;*/
	double tau = 0.05958;


	// solver
	TestTransport3WorkerExplicit worker(n_cells);

	//V = worker.V0();
	//U = worker.U0();

	// saver
	VtkUtils::TimeDependentWriter writer("transport3-explicit");
	std::string out_filename = writer.add(0);
	worker.save_vtk(out_filename);

	double norm = 0;

	while (worker.current_time() < tend - 1e-6) {
		// solve problem
		norm = worker.step(tau);
		// export solution to vtk
		out_filename = writer.add(worker.current_time());
		worker.save_vtk(out_filename);
	};
	std::cout << 1.0 / tau << " " << norm << std::endl;
	//CHECK(norm == Approx(0.0138123932).margin(1e-5));
}