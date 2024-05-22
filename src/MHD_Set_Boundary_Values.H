#pragma once

#include "Proto.H"
#include "MHDOp.H"
#include "MHDLevelDataRK4.H"
namespace MHD_Set_Boundary_Values {

	void Set_Boundary_Values_Spherical_2O(LevelBoxData<double,NUMCOMPS>& a_U,
	                         MHDLevelDataState& a_State);					 

	void interpolate_h5_BC(MHDLevelDataState& state,
						const BoxData<double,NUMCOMPS>& BC_data,
						const double time);

	void interpolate_h5_BC(MHDLevelDataState& state,
						const std::vector<BoxData<double,NUMCOMPS>>& BC_data,
						const double time);					
}
