#pragma once

#include "Proto.H"
#include "MHDOp.H"

using namespace Proto;

namespace MHD_Initialize
{
	void initializeState(LevelBoxData<double, NUMCOMPS> &a_U,
						 const double dx,
						 const double dy,
						 const double dz);

	void initializeState(AMRData<double, NUMCOMPS> &a_U,
						 const double dx,
						 const double dy,
						 const double dz);

	void initializeState(MHDLevelDataState &a_State);

	void initializeState_Spherical(MHDLevelDataState &a_State);

	void InitializeStatecalc(BoxData<double, NUMCOMPS> &a_U,
							 const BoxData<double, DIM> &a_x,
							 const BoxData<double, DIM> &a_eta,
							 const double a_gamma,
							 const Box &a_dbx);

	void InitializeState_Spherical_2O(BoxData<double, NUMCOMPS> &a_U,
									  const BoxData<double, DIM> &a_x,
									  const double a_gamma);

	void initializeState_Spherical_2O(MHDLevelDataState &a_State);

}