set (SRC
	cfd24.cpp

	mat/csrmat.cpp
	mat/lodmat.cpp
	mat/sparse_matrix_solver.cpp

	grid/vtk.cpp
	grid/grid1d.cpp
	grid/regular_grid2d.cpp

	debug/printer.cpp
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

	grid/vtk.hpp
	grid/i_grid.hpp
	grid/grid1d.hpp
	grid/regular_grid2d.hpp

	geom/primitives.hpp

	debug/printer.hpp
)

add_library(${CFD24_TARGET} STATIC ${SRC} ${HEADERS})
