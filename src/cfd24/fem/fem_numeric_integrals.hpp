#ifndef __CFD_FEM_NUMERIC_INTEGRALS_HPP__
#define __CFD_FEM_NUMERIC_INTEGRALS_HPP__

#include "cfd24/fem/fem_element.hpp"
#include "cfd24/numeric_integration/quadrature.hpp"

namespace cfd{

class NumericElementIntegrals: public IElementIntegrals{
public:
	NumericElementIntegrals(
			const Quadrature* quad,
			std::shared_ptr<const IElementGeometry> geom,
			std::shared_ptr<const IElementBasis> basis);

	std::vector<double> mass_matrix() const override;
	std::vector<double> load_vector() const override;
	std::vector<double> stiff_matrix() const override;
	std::vector<double> transport_matrix(
			const std::vector<double>& vx,
			const std::vector<double>& vy={},
			const std::vector<double>& vz={}) const override;

	std::vector<double> mass_matrix_stab_supg(
			const std::vector<double>& vx,
			const std::vector<double>& vy={},
			const std::vector<double>& vz={}) const override;
	std::vector<double> load_vector_stab_supg(
			const std::vector<double>& vx,
			const std::vector<double>& vy={},
			const std::vector<double>& vz={}) const override;
	std::vector<double> stiff_matrix_stab_supg(
			const std::vector<double>& vx,
			const std::vector<double>& vy={},
			const std::vector<double>& vz={}) const override;
	std::vector<double> transport_matrix_stab_supg(
			const std::vector<double>& vx,
			const std::vector<double>& vy={},
			const std::vector<double>& vz={}) const override;
private:
	const Quadrature* _quad;
	std::shared_ptr<const IElementGeometry> _geom;
	std::shared_ptr<const IElementBasis> _basis;
	const std::vector<JacobiMatrix> _quad_jacobi;
};


}
#endif
