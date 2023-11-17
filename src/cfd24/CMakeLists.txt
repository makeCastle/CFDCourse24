set (SRC
	cfd24.cpp

	mat/csrmat.cpp
	mat/lodmat.cpp
	mat/densemat.cpp
	mat/sparse_matrix_solver.cpp

	geom/simplex.cpp
	geom/searcher.cpp

	grid/vtk.cpp
	grid/grid1d.cpp
	grid/regular_grid2d.cpp
	grid/unstructured_grid2d.cpp

	fvm/fvm_assembler.cpp

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

	geom/primitives.hpp
	geom/simplex.hpp
	geom/searcher.hpp

	grid/vtk.hpp
	grid/i_grid.hpp
	grid/grid1d.hpp
	grid/regular_grid2d.hpp
	grid/unstructured_grid2d.hpp

	fvm/fvm_assembler.hpp

	debug/printer.hpp
	debug/saver.hpp
	debug/tictoc.hpp
)

add_library(${CFD24_TARGET} STATIC ${SRC} ${HEADERS})
