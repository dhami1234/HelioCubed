#pragma once

#include "Proto.H"
#include "MHDOp.H"
//TODO: These variables need to go into a common file or be defined in a more accessible way
//constexpr unsigned int NUMCOMPS=2*DIM+2;
typedef Var<double,NUMCOMPS> State;


namespace MHD_Limiters {

	void MHD_Limiters_minmod(BoxData<double,NUMCOMPS>& a_W_low_lim,
	                  BoxData<double,NUMCOMPS>& a_W_high_lim,
	                  BoxData<double,NUMCOMPS>& a_W,
	                  BoxData<double,DIM>& a_x_sph,
	                  BoxData<double,DIM>& a_dx_sph,
	                  const int a_d);				  
}
