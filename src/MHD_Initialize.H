#pragma once

#include "Proto.H"
#include "MHDOp.H"

using namespace Proto;

namespace MHD_Initialize
{
					 
	void InitializeState_Spherical_2O(BoxData<double,NUMCOMPS>& a_U,
	                         const BoxData<double,DIM>& a_x,
	                         const double a_gamma);


	void initializeState_Spherical_2O(MHDLevelDataState& a_State);
						 				 		 
}