set (SRC
	cfd24.cpp

	grid/vtk.cpp
	grid/grid1d.cpp
)

set (HEADERS
	cfd24.hpp

	grid/vtk.hpp
	grid/i_grid.hpp
	grid/grid1d.hpp

	geom/primitives.hpp
)

add_library(${CFD24_TARGET} STATIC ${SRC} ${HEADERS})
