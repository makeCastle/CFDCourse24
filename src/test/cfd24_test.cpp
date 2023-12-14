#define CATCH_CONFIG_RUNNER
#include "cfd24_test.hpp"
#include "cfd24/debug/saver.hpp"
#include "cfd24/debug/printer.hpp"

TEST_CASE("Ping", "[ping]"){
	cfd::dbg::ping_saver_cpp();
	cfd::dbg::ping_printer_cpp();
	CHECK(cfd::ping() == 1);
}


int main(int argc, char* argv[]){
	int result = Catch::Session().run(argc, argv);
	std::cout << "DONE" << std::endl;
	return result;
}
