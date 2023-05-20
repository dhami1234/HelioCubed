#pragma once

#include "Proto.H"

#include "MHDOp.H"

using namespace Proto;

namespace MHD_Artificial_Viscosity
{

	// Used to implement artificial viscosity
	void step(LevelBoxData<double, NUMCOMPS> &a_Rhs,
			  LevelBoxData<double, NUMCOMPS> &a_JU,
			  MHDLevelDataState &a_State);

	void step_spherical(LevelBoxData<double, NUMCOMPS> &a_Rhs,
						LevelBoxData<double, NUMCOMPS> &a_JU,
						MHDLevelDataState &a_State);
}
