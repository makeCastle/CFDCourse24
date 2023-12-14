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

const Quadrature* cfd::quadrature_triangle_gauss5(){
	static Quadrature quad(
		{
			Point(0.33333333333333, 0.33333333333333),
			Point(0.47014206410511, 0.47014206410511),
			Point(0.47014206410511, 0.05971587178977),
			Point(0.05971587178977, 0.47014206410511),
			Point(0.10128650732346, 0.10128650732346),
			Point(0.10128650732346, 0.79742698535309),
			Point(0.79742698535309, 0.10128650732346)
		},
		{
			0.22500000000000/2,
			0.13239415278851/2,
			0.13239415278851/2,
			0.13239415278851/2,
			0.12593918054483/2,
			0.12593918054483/2,
			0.12593918054483/2
		});
	return &quad;
}

const Quadrature* cfd::quadrature_triangle_gauss6(){
	static Quadrature quad(
		{
			Point(0.24928674517091, 0.24928674517091),
			Point(0.24928674517091, 0.50142650965818),
			Point(0.50142650965818, 0.24928674517091),
			Point(0.06308901449150, 0.06308901449150),
			Point(0.06308901449150, 0.87382197101700),
			Point(0.87382197101700, 0.06308901449150),
			Point(0.31035245103378, 0.63650249912140),
			Point(0.63650249912140, 0.05314504984482),
			Point(0.05314504984482, 0.31035245103378),
			Point(0.63650249912140, 0.31035245103378),
			Point(0.31035245103378, 0.05314504984482),
			Point(0.05314504984482, 0.63650249912140),
		},
		{
			0.11678627572638/2,
			0.11678627572638/2,
			0.11678627572638/2,
			0.05084490637021/2,
			0.05084490637021/2,
			0.05084490637021/2,
			0.08285107561837/2,
			0.08285107561837/2,
			0.08285107561837/2,
			0.08285107561837/2,
			0.08285107561837/2,
			0.08285107561837/2
		});
	return &quad;
}
