#include "jacobi.hpp"

using namespace cfd;

void cfd::fill_jacobi_modj(JacobiMatrix& jac){
	jac.modj = jac.j11*(jac.j22*jac.j33 - jac.j23*jac.j32)
	          -jac.j12*(jac.j21*jac.j33 - jac.j23*jac.j31)
	          +jac.j13*(jac.j21*jac.j32 - jac.j22*jac.j31);
}

void cfd::fill_jacobi_modj_1d(JacobiMatrix& jac){
	jac.modj = jac.j11;
}

void cfd::fill_jacobi_modj_2d(JacobiMatrix& jac){
	jac.modj = jac.j11 * jac.j22 - jac.j12*jac.j21;
}

Vector cfd::gradient_to_parametric(const JacobiMatrix& jac, Vector grad_x){
	return {
		jac.j11 * grad_x.x() + jac.j21 * grad_x.y() + jac.j31 * grad_x.z(),
		jac.j12 * grad_x.x() + jac.j22 * grad_x.y() + jac.j32 * grad_x.z(),
		jac.j13 * grad_x.x() + jac.j23 * grad_x.y() + jac.j33 * grad_x.z(),
	};
}

Vector cfd::gradient_to_parametric_1d(const JacobiMatrix& jac, Vector grad_x){
	return {
		jac.j11 * grad_x.x()
	};
}

Vector cfd::gradient_to_parametric_2d(const JacobiMatrix& jac, Vector grad_x){
	return {
		jac.j11 * grad_x.x() + jac.j21 * grad_x.y(),
		jac.j12 * grad_x.x() + jac.j22 * grad_x.y()
	};
}

Vector cfd::gradient_to_physical(const JacobiMatrix& jac, Vector grad_xi){
	double c = 1.0/jac.modj;
	double a11 = c*(jac.j22*jac.j33-jac.j23*jac.j32);
	double a12 = c*(jac.j23*jac.j31-jac.j21*jac.j33);
	double a13 = c*(jac.j21*jac.j32-jac.j22*jac.j31);
	double a21 = c*(jac.j13*jac.j32-jac.j12*jac.j33);
	double a22 = c*(jac.j11*jac.j33-jac.j13*jac.j31);
	double a23 = c*(jac.j12*jac.j31-jac.j11*jac.j32);
	double a31 = c*(jac.j12*jac.j23-jac.j13*jac.j22);
	double a32 = c*(jac.j13*jac.j21-jac.j11*jac.j23);
	double a33 = c*(jac.j11*jac.j22-jac.j12*jac.j21);

	return {
		a11 * grad_xi.x() + a12 * grad_xi.y() + a13 * grad_xi.z(),
		a21 * grad_xi.x() + a22 * grad_xi.y() + a23 * grad_xi.z(),
		a31 * grad_xi.x() + a32 * grad_xi.y() + a33 * grad_xi.z()
	};
}

Vector cfd::gradient_to_physical_1d(const JacobiMatrix& jac, Vector grad_xi){
	return { grad_xi.x() / jac.j11 };
}

Vector cfd::gradient_to_physical_2d(const JacobiMatrix& jac, Vector grad_xi){
	double c = 1.0/jac.modj;
	double a11 = c*jac.j22;
	double a12 = -c*jac.j21;
	double a21 = -c*jac.j12;
	double a22 = c*jac.j11;

	return {
		a11 * grad_xi.x() + a12 * grad_xi.y(),
		a21 * grad_xi.x() + a22 * grad_xi.y()
	};
}

double cfd::laplace_to_physical(const JacobiMatrix& jac, const std::array<double, 6>& upper_hessian_xi){
	double c = 1.0/jac.modj;
	double a11 = c*(jac.j22*jac.j33-jac.j23*jac.j32);
	double a12 = c*(jac.j23*jac.j31-jac.j21*jac.j33);
	double a13 = c*(jac.j21*jac.j32-jac.j22*jac.j31);
	double a21 = c*(jac.j13*jac.j32-jac.j12*jac.j33);
	double a22 = c*(jac.j11*jac.j33-jac.j13*jac.j31);
	double a23 = c*(jac.j12*jac.j31-jac.j11*jac.j32);
	double a31 = c*(jac.j12*jac.j23-jac.j13*jac.j22);
	double a32 = c*(jac.j13*jac.j21-jac.j11*jac.j23);
	double a33 = c*(jac.j11*jac.j22-jac.j12*jac.j21);

	const auto& h = upper_hessian_xi;

	double dx2 = a11*a11*h[0] + a12*a12*h[3] + a13*a13*h[5] + 2*a11*a12*h[1] + 2*a11*a13*h[2] + 2*a12*a13*h[4];
	double dy2 = a21*a21*h[0] + a22*a22*h[3] + a23*a23*h[5] + 2*a21*a22*h[1] + 2*a21*a23*h[2] + 2*a22*a23*h[4];
	double dz2 = a31*a31*h[0] + a32*a32*h[3] + a33*a33*h[5] + 2*a31*a32*h[1] + 2*a31*a33*h[2] + 2*a32*a33*h[4];
	
	return dx2 + dy2 + dz2;
}
