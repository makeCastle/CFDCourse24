#ifndef UTILS_FILESYSTEM_HPP
#define UTILS_FILESYSTEM_HPP

#include <string>

/**
 * Creates path to the file located in the test_data directory
 *
 * @param path  relative path to the file in the test data directory
 *
 * @return      full path
 */
std::string test_directory(std::string path);


#endif