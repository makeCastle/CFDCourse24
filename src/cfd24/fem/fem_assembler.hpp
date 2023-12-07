#ifndef __CFD_FEM_ASSEMBLER_HPP__
#define __CFD_FEM_ASSEMBLER_HPP__

#include "cfd24/fem/fem_element.hpp"
#include "cfd24/mat/csrmat.hpp"
#include "cfd24/grid/i_grid.hpp"

namespace cfd{

struct FemAssembler{
	FemAssembler(size_t n_bases,
	             const std::vector<FemElement>& elements,
	             const std::vector<std::vector<size_t>>& tab_elem_basis);

	virtual ~FemAssembler() = default;

	size_t n_elements() const;
	size_t n_bases() const;
	const CsrStencil& stencil() const;

	const FemElement& element(size_t icell) const;

	Point reference_point(size_t ibas) const;

	void add_to_global_matrix(size_t ielem, const std::vector<double>& local_matrix, std::vector<double>& global_csr_vals) const;
	void add_to_global_vector(size_t ielem, const std::vector<double>& local_vector, std::vector<double>& global_vector) const;
protected:
	std::vector<FemElement> _elements;

	std::vector<std::vector<size_t>> _tab_elem_basis;
	std::vector<std::vector<size_t>> _tab_elem_csr_address;
	std::vector<Point> _ref_points;

	CsrStencil _stencil;
};

}
#endif
