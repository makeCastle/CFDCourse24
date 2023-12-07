set (SRC
	cfd24.cpp

	mat/csrmat.cpp
	mat/lodmat.cpp
	mat/densemat.cpp
	mat/sparse_matrix_solver.cpp

	geom/simplex.cpp
	geom/searcher.cpp
	geom/jacobi.cpp

	grid/vtk.cpp
	grid/i_grid.cpp
	grid/grid1d.cpp
	grid/regular_grid2d.cpp
	grid/unstructured_grid2d.cpp

	numeric_integration/quadrature.cpp
	numeric_integration/segment_quadrature.cpp
	numeric_integration/square_quadrature.cpp
	numeric_integration/triangle_quadrature.cpp

	fvm/fvm_assembler.cpp
	fvm/fvm_dpdn_boundary.cpp

	fem/fem_assembler.cpp
	fem/fem_numeric_integrals.cpp
	fem/elem1d/segment_linear.cpp
	fem/elem2d/triangle_linear.cpp
	fem/elem2d/quadrangle_linear.cpp

	debug/printer.cpp
	debug/saver.cpp
	debug/tictoc.cpp
)

set (HEADERS
	cfd24.hpp
	cfd_common.hpp
	macros.hpp

	mat/i_mat.hpp
	mat/i_sparse_mat.hpp
	mat/csrmat.hpp
	mat/lodmat.hpp
	mat/sparse_matrix_solver.hpp
	mat/densemat.hpp

	geom/point.hpp
	geom/simplex.hpp
	geom/searcher.hpp
	geom/jacobi.hpp

	grid/vtk.hpp
	grid/i_grid.hpp
	grid/grid1d.hpp
	grid/regular_grid2d.hpp
	grid/unstructured_grid2d.hpp

	numeric_integration/quadrature.hpp
	numeric_integration/segment_quadrature.hpp
	numeric_integration/square_quadrature.hpp
	numeric_integration/triangle_quadrature.hpp

	fvm/fvm_assembler.hpp
	fvm/fvm_dpdn_boundary.hpp

	fem/fem_assembler.hpp
	fem/fem_element.hpp
	fem/fem_numeric_integrals.hpp
	fem/elem1d/segment_linear.hpp
	fem/elem2d/triangle_linear.hpp
	fem/elem2d/quadrangle_linear.hpp

	debug/printer.hpp
	debug/saver.hpp
	debug/tictoc.hpp
)

add_library(${CFD24_TARGET} STATIC ${SRC} ${HEADERS})
