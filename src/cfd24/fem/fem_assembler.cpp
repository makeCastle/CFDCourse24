#include "fem_assembler.hpp"

using namespace cfd;

FemAssembler::FemAssembler(size_t n_bases, const std::vector<FemElement>& elements,
		const std::vector<std::vector<size_t>>& tab_elem_basis):
		_elements(elements), _tab_elem_basis(tab_elem_basis){
	
	// stencil
	std::vector<std::set<size_t>> tab_basis_basis(n_bases);
	for (size_t ielem=0; ielem<elements.size(); ++ielem)
	for (size_t i=0; i<tab_elem_basis[ielem].size(); ++i)
	for (size_t j=0; j<i+1; ++j){
		size_t ibas1 = tab_elem_basis[ielem][i];
		size_t ibas2 = tab_elem_basis[ielem][j];
		tab_basis_basis[ibas1].insert(ibas2);
		tab_basis_basis[ibas2].insert(ibas1);
	}
	_stencil.set_stencil(tab_basis_basis);

	// element -> csr addresses table
	_tab_elem_csr_address.resize(n_elements());
	for (size_t ielem=0; ielem<n_elements(); ++ielem){
		size_t nbas = _elements[ielem].basis->size();
		for (size_t local_row=0; local_row<nbas; ++local_row){
			size_t global_row = tab_elem_basis[ielem][local_row];
			for (size_t local_col=0; local_col<nbas; ++local_col){
				size_t global_col = tab_elem_basis[ielem][local_col];
				size_t addr = _stencil.get_address(global_row, global_col);
				if (addr == INVALID_INDEX){
					_THROW_INTERNAL_ERROR_;
				}
				_tab_elem_csr_address[ielem].push_back(addr);
			}
		}
	}

	// basis types and reference points
	_ref_points.resize(n_bases);
	_bas_types.resize(n_bases);
	for (size_t ielem=0; ielem < elements.size(); ++ielem){
		auto geom = _elements[ielem].geometry;
		std::vector<Point> xi = _elements[ielem].basis->parametric_reference_points();
		std::vector<BasisType> types = _elements[ielem].basis->basis_types();
		for (size_t ibas = 0; ibas < xi.size(); ++ibas){
			Point p = geom->to_physical(xi[ibas]);
			size_t ind = _tab_elem_basis[ielem][ibas];
			_ref_points[ind] = p;
			_bas_types[ind] = types[ibas];
		}
	}
}

size_t FemAssembler::n_elements() const{
	return _elements.size();
}

size_t FemAssembler::n_bases() const{
	return _stencil.n_rows();
}

const FemElement& FemAssembler::element(size_t ielem) const{
	return _elements[ielem];
}

Point FemAssembler::reference_point(size_t ibas) const{
	return _ref_points[ibas];
}

std::vector<double> FemAssembler::approximate(const IPointFunction& func) const{
	std::vector<double> ret;
	for (size_t ibas=0; ibas < n_bases(); ++ibas){
		Point p = reference_point(ibas);
		switch (_bas_types[ibas]){
			case BasisType::Nodal: ret.push_back(func(p)); break;
			case BasisType::Dx: ret.push_back(func.grad(p).x()); break;
			case BasisType::Dy: ret.push_back(func.grad(p).y()); break;
			case BasisType::Dz: ret.push_back(func.grad(p).z()); break;
			default: _THROW_NOT_IMP_;
		};
	}
	return ret;
}

std::vector<double> FemAssembler::local_vector(size_t ielem, const std::vector<double>& v) const{
	std::vector<double> ret;
	for (size_t bas: _tab_elem_basis[ielem]){
		ret.push_back(v[bas]);
	}
	return ret;
}

void FemAssembler::add_to_global_matrix(size_t ielem, const std::vector<double>& local_matrix, std::vector<double>& global_csr_vals) const{
	for (size_t ival=0; ival < local_matrix.size(); ++ival){
		size_t a = _tab_elem_csr_address[ielem][ival];
		global_csr_vals[a] += local_matrix[ival];
	}
}

void FemAssembler::add_to_global_vector(size_t ielem, const std::vector<double>& local_vector, std::vector<double>& global_vector) const{
	for (size_t i=0; i<local_vector.size(); ++i){
		size_t gi = _tab_elem_basis[ielem][i];
		global_vector[gi] += local_vector[i];
	}
}

const CsrStencil& FemAssembler::stencil() const{
	return _stencil;
}
