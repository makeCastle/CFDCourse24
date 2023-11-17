#include "fvm_assembler.hpp"
#include "cfd24/geom/searcher.hpp"
#include "cfd24/geom/simplex.hpp"
#include "cfd24/debug/tictoc.hpp"
#include "cfd24/mat/densemat.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// FvmExtendedCollocations
///////////////////////////////////////////////////////////////////////////////

FvmExtendedCollocations::FvmExtendedCollocations(const IGrid& grid){
	// cell loop
	for (size_t icell=0; icell<grid.n_cells(); ++icell){
		// collocations at the cell center
		points.push_back(grid.cell_center(icell));
		cell_collocations.push_back({points.size()-1, icell});
	}

	// face loop
	for (size_t iface=0; iface<grid.n_faces(); ++iface){
		std::array<size_t, 2> cells = grid.tab_face_cell(iface);

		// face -> collocation-points connectivity
		tab_face_colloc.push_back({cells[0], cells[1]});

		// collocations at the boundary faces
		if (cells[0] == INVALID_INDEX || cells[1] == INVALID_INDEX){
			points.push_back(grid.face_center(iface));
			face_collocations.push_back({points.size()-1, iface});

			if (cells[0] == INVALID_INDEX){
				tab_face_colloc.back().negative_side = points.size()-1;
			} else {
				tab_face_colloc.back().positive_side = points.size()-1;
			}
		}
	}

	// collocations connectivity
	tab_colloc_colloc.resize(points.size());
	for (const FaceConnect& fc: tab_face_colloc){
		tab_colloc_colloc[fc.positive_side].push_back(fc.negative_side);
		tab_colloc_colloc[fc.negative_side].push_back(fc.positive_side);
	}
}

size_t FvmExtendedCollocations::size() const{
	return points.size();
}

///////////////////////////////////////////////////////////////////////////////
// Cell gradient
///////////////////////////////////////////////////////////////////////////////

namespace{

DenseMatrix least_squares_inv(const DenseMatrix& a){
	// transpose(A)
	DenseMatrix at = a.transpose();
	// inverse(transpose(A) * A)
	DenseMatrix inv_at_a = at.mult_mat(a).inverse();
	// inverse(transpose(A) * A) * transpose(A)
	return inv_at_a.mult_mat(at);
}

std::array<CsrMatrix, 3> assemble_fvm_cell_gradient_2d(const IGrid& grid, const FvmExtendedCollocations& colloc){
	LodMatrix grad_x(grid.n_cells());
	LodMatrix grad_y(grid.n_cells());

	for (size_t icell = 0; icell < grid.n_cells(); ++icell){
		const std::vector<size_t>& collocs = colloc.tab_colloc_colloc[icell];

		DenseMatrix amat(collocs.size(), 2);
		for (size_t i=0; i<collocs.size(); ++i){
			Vector c = colloc.points[collocs[i]] - colloc.points[icell];
			amat.set_value(i, 0, c.x());
			amat.set_value(i, 1, c.y());
		}
		DenseMatrix lsi = least_squares_inv(amat);
		double diag_x = 0;
		double diag_y = 0;
		for (size_t i=0; i<collocs.size(); ++i){
			double vx = lsi.value(0, i);
			double vy = lsi.value(1, i);
			grad_x.set_value(icell, collocs[i], vx);
			grad_y.set_value(icell, collocs[i], vy);
			diag_x -= vx;
			diag_y -= vy;
		}
		grad_x.set_value(icell, icell, diag_x);
		grad_y.set_value(icell, icell, diag_y);
	}

	return {grad_x.to_csr(), grad_y.to_csr()};
}

std::array<CsrMatrix, 3> assemble_fvm_cell_gradient(const IGrid& grid, const FvmExtendedCollocations& colloc){
	if (grid.dim() == 2){
		return assemble_fvm_cell_gradient_2d(grid, colloc);
	} else {
		_THROW_NOT_IMP_;
	}
}

}

FvmCellGradient::FvmCellGradient(const IGrid& grid, const FvmExtendedCollocations& colloc)
	: _data(assemble_fvm_cell_gradient(grid, colloc)){}

std::vector<Vector> FvmCellGradient::compute(const std::vector<double>& u) const{
	std::vector<double> x = _data[0].mult_vec(u);
	std::vector<double> y = _data[1].mult_vec(u);

	std::vector<Vector> ret(x.size());
	for (size_t i=0; i<ret.size(); ++i){
		ret[i].x() = x[i];
		ret[i].y() = y[i];
	}

	if (_data[2].n_rows() > 0){
		std::vector<double> z = _data[2].mult_vec(u);
		for (size_t i=0; i<ret.size(); ++i){
			ret[i].z() = z[i];
		}
	}

	return ret;
}

///////////////////////////////////////////////////////////////////////////////
// DuDn
///////////////////////////////////////////////////////////////////////////////

namespace {

LodMatrix assemble_faces_dudn_2d(const IGrid& grid, const FvmExtendedCollocations& colloc){
	LodMatrix mat(grid.n_faces());

	std::vector<std::array<size_t, 3>> closest_points;
	{
		PointSearcher<2> searcher(colloc.points);
		for (size_t ipoint=0; ipoint < grid.n_points(); ++ipoint){
			std::vector<size_t> found = searcher.nearest(grid.point(ipoint), 3);
			closest_points.push_back({found[0], found[1], found[2]});
		}
	}
	auto find_closest_collocation = [&closest_points](size_t grid_point, size_t excl0, size_t excl1){
		const std::array<size_t, 3>& found = closest_points[grid_point];
		if (found[0] != excl0 && found[0] != excl1){
			return found[0];
		} else if (found[1] != excl0 && found[1] != excl1){
			return found[1];
		} else {
			return found[2];
		}
	};

	auto add_ds_entry = [&](const Vector& normal, const Vector& c, size_t col0, size_t col1, size_t col2, size_t iface){
		Vector s(-normal.y(), normal.x());
		double cos_cos = dot_product(c, s) / dot_product(c, normal);

		Point p0 = colloc.points[col0];
		Point p1 = colloc.points[col1];
		Point p2 = colloc.points[col2];
		double tri_area = triangle_area(p0, p1, p2);

		double coef = -0.5*cos_cos/tri_area;

		double x0 = p0.x(); double y0 = p0.y();
		double x1 = p1.x(); double y1 = p1.y();
		double x2 = p2.x(); double y2 = p2.y();
		double dx0 = (y1 - y2)/2.0;
		double dy0 = (x2 - x1)/2.0;
		double dx1 = (y2 - y0)/2.0;
		double dy1 = (x0 - x2)/2.0;
		double dx2 = (y0 - y1)/2.0;
		double dy2 = (x1 - x0)/2.0;

		mat.add_value(iface, col0, coef*(dx0*s.x() + dy0*s.y()));
		mat.add_value(iface, col1, coef*(dx1*s.x() + dy1*s.y()));
		mat.add_value(iface, col2, coef*(dx2*s.x() + dy2*s.y()));
	};

	for (size_t iface = 0; iface < grid.n_faces(); ++iface){
		Vector normal = grid.face_normal(iface);
		size_t negative_collocation = colloc.tab_face_colloc[iface].negative_side;
		size_t positive_collocation = colloc.tab_face_colloc[iface].positive_side;
		Point ci = colloc.points[negative_collocation];
		Point cj = colloc.points[positive_collocation];
	
		// +dudc / cos(c,n);
		double v1 = 1.0/dot_product(normal, cj-ci);
		mat.set_value(iface, positive_collocation, v1);
		mat.set_value(iface, negative_collocation, -v1);

		// -duds*cos(c,s)/cos(c,n) 
		Vector c = (cj - ci)/vector_abs(cj - ci);
		{
			// left point
			size_t igrid = grid.tab_face_point(iface)[0];
			size_t col1 = find_closest_collocation(igrid, positive_collocation, negative_collocation);
			add_ds_entry(normal, c, negative_collocation, col1, positive_collocation, iface);
		}
		{
			// right point
			size_t igrid = grid.tab_face_point(iface)[1];
			size_t col1 = find_closest_collocation(igrid, positive_collocation, negative_collocation);
			add_ds_entry(normal, c, positive_collocation, col1, negative_collocation, iface);
		}
	}

	return mat;
}

}

LodMatrix cfd::assemble_fvm_faces_dudn(const IGrid& grid, const FvmExtendedCollocations& colloc){
	if (grid.dim() == 2){
		return assemble_faces_dudn_2d(grid, colloc);
	} else {
		_THROW_NOT_IMP_;
	}
}

