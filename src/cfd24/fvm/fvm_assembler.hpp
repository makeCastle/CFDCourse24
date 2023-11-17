#ifndef CFD_FVM_ASSEMBLER_HPP
#define CFD_FVM_ASSEMBLER_HPP

#include "cfd24/mat/lodmat.hpp"
#include "cfd24/grid/i_grid.hpp"

namespace cfd{

///////////////////////////////////////////////////////////////////////////////
// FvmExtendedCollocations
///////////////////////////////////////////////////////////////////////////////

struct FvmExtendedCollocations{
	FvmExtendedCollocations(const IGrid& grid);

	std::vector<Point> points;

	size_t size() const;

	struct FaceConnect{
		size_t negative_side;
		size_t positive_side;
	};
	std::vector<FaceConnect> tab_face_colloc;

	struct IndexConnect{
		size_t colloc_index;
		size_t grid_index;
	};
	std::vector<IndexConnect> cell_collocations;
	std::vector<IndexConnect> face_collocations;

	std::vector<std::vector<size_t>> tab_colloc_colloc;
};

///////////////////////////////////////////////////////////////////////////////
// FvmCellGradient
///////////////////////////////////////////////////////////////////////////////

struct FvmCellGradient{
	FvmCellGradient(const IGrid& grid, const FvmExtendedCollocations& colloc);

	std::vector<Vector> compute(const std::vector<double>& u) const;
private:
	std::array<CsrMatrix, 3> _data;
};

///////////////////////////////////////////////////////////////////////////////
// DuDn on faces
///////////////////////////////////////////////////////////////////////////////

LodMatrix assemble_fvm_faces_dudn(const IGrid& grid, const FvmExtendedCollocations& colloc);

}
#endif
