#include "utils/filesystem.hpp"
#include <iostream>

#ifndef TEST_DIRECTORY
#define TEST_DIRECTORY "./"
#endif

std::string test_directory(std::string path){
	return TEST_DIRECTORY + path;
}

