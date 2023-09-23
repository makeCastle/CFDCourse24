#define CATCH_CONFIG_RUNNER
#include "cfd24_test.hpp"

#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/grid/vtk.hpp"

double func(double x, double y, double t){
	//x = x - t;
	//if (x >= -1 && x <=-0.8 && y >= -0.1 && y <= 0.1){
	//        return 1;
	//} else {
	//        return 0;
	//}
	double sigma = 0.1;
	
	//double vx = -y;
	//double vy = x;

	double x0 = 0.5;
	double y0 = 0;

	double xc =  x0*cos(t/2);
	double yc =  x0*sin(t/2);

	double r = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));

	return exp(-r*r/sigma/sigma);
}

TEST_CASE("Ping", "[ping]"){
	CHECK(cfd::ping() == 1);
	cfd::RegularGrid2D grid(-1, 1, -1, 1, 100, 100);

	std::vector<double> v(grid.n_points());
	for (size_t i=0; i<grid.n_points(); ++i){
		v[i] = func(grid.point(i).x(), grid.point(i).y(), M_PI);
	}

	grid.save_vtk("out.vtk");
	cfd::VtkUtils::add_point_data(v, "data", "out.vtk");
}

int main(int argc, char* argv[]){
	int result = Catch::Session().run(argc, argv);
	std::cout << "DONE" << std::endl;
	return result;
}
