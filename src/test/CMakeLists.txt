set (SRC
	cfd24_test.cpp
	grid_test.cpp
	mat_test.cpp
	poisson_solve_test.cpp
	transport_solve_test.cpp
	cavern_2d_simple_test.cpp

	utils/filesystem.cpp
)

set (HEADERS
	cfd24_test.hpp

	utils/filesystem.hpp
)

add_executable(cfd24_test ${SRC} ${HEADERS})
target_link_libraries(cfd24_test PRIVATE ${CATCH2_FRAMEWORK} ${CFD24_TARGET})
target_include_directories(cfd24_test PRIVATE ${CMAKE_SOURCE_DIR}/src/test)
target_compile_definitions(cfd24_test PUBLIC TEST_DIRECTORY="${CMAKE_SOURCE_DIR}/test_data/")
