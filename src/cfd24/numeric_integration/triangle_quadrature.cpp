#include "triangle_quadrature.hpp"

using namespace cfd;

const Quadrature* cfd::quadrature_triangle_gauss1(){
	static Quadrature quad(
		{ Point(1.0/3.0, 1.0/3.0) },
		{ 0.5 });

	return &quad;
}

const Quadrature* cfd::quadrature_triangle_gauss2(){
	static Quadrature quad(
		{ 
			Point(1.0/6.0, 1.0/6.0),
			Point(2.0/3.0, 1.0/6.0),
			Point(1.0/6.0, 2.0/3.0),
		},
		{ 
			1.0/6.0,
			1.0/6.0,
			1.0/6.0
		});

	return &quad;
}

const Quadrature* cfd::quadrature_triangle_gauss3(){
	static Quadrature quad(
		{ 
			Point(1.0/3.0, 1.0/3.0),
			Point(1.0/5.0, 1.0/5.0),
			Point(1.0/5.0, 3.0/5.0),
			Point(3.0/5.0, 1.0/5.0),
		},
		{ 
			-27.0/96.0,
			 25.0/96.0,
			 25.0/96.0,
			 25.0/96.0
		});

	return &quad;
}

const Quadrature* cfd::quadrature_triangle_gauss4(){
	static Quadrature quad(
		{ 
			Point(0.44594849091597, 0.44594849091597),
			Point(0.44594849091597, 0.10810301816807),
			Point(0.10810301816807, 0.44594849091597),
			Point(0.09157621350977, 0.09157621350977),
			Point(0.09157621350977, 0.81684757298046),
			Point(0.81684757298046, 0.09157621350977)
		},
		{ 
			0.22338158967801/2.0,
			0.22338158967801/2.0,
			0.22338158967801/2.0,
			0.10995174365532/2.0,
			0.10995174365532/2.0,
			0.10995174365532/2.0
		});

	return &quad;
}
