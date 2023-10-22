#pragma once

#include "Proto.H"

#include "MHDOp.H"

using namespace Proto;

namespace MHD_Turbulence {


	void Turb_Flux(BoxData<double,NUMCOMPS>& a_F,
			  BoxData<double,NUMCOMPS>& a_W_lo,
              BoxData<double,NUMCOMPS>& a_W_hi,
			  const int a_dir);

	void Turb_Source(BoxData<double,NUMCOMPS>& a_S,
			  BoxData<double,NUMCOMPS>& a_W,
              BoxData<double,1>& a_divV);

}