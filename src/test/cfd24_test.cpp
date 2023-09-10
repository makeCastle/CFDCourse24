#define CATCH_CONFIG_RUNNER
#include "cfd24_test.hpp"

#include "cfd24/grid/regular_grid2d.hpp"
#include "cfd24/grid/vtk.hpp"

TEST_CASE("Ping", "[ping]"){
	CHECK(cfd::ping() == 1);
}

int main(int argc, char* argv[]){
	int result = Catch::Session().run(argc, argv);
	std::cout << "DONE" << std::endl;
	return result;
}
