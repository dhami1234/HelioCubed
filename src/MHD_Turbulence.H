#pragma once

#include "Proto.H"

#include "MHDOp.H"

using namespace Proto;

namespace MHD_Turbulence {

	/// Physical flux used by the MUSCL-Hancock predictor. As in MS-FLUKSS,
	/// this only advances the three turbulence conserved variables.
	void Turb_Flux_Hancock(BoxData<double,NUMCOMPS>& a_F,
			  BoxData<double,NUMCOMPS>& a_W,
			  const int a_dir);

	void Turb_Flux(BoxData<double,NUMCOMPS>& a_F,
			  BoxData<double,NUMCOMPS>& a_W_lo,
              BoxData<double,NUMCOMPS>& a_W_hi,
			  const int a_dir);

	void Turb_Source(BoxData<double,NUMCOMPS>& a_S,
			  BoxData<double,NUMCOMPS>& a_W,
              BoxData<double,1>& a_divV);

}
