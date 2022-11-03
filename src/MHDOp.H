#pragma once

#include "Proto.H"

//TODO: These variables need to go into a common file or be defined in a more accessible way
constexpr unsigned int NUMCOMPS=2*DIM+2;

#define NGHOST 2

using namespace Proto;

typedef Var<double,DIM> V;
typedef Var<double,NUMCOMPS> State;
class MHDLevelDataState;

namespace MHDOp
{
	void DimToNonDimcalc(BoxData<double,NUMCOMPS>& a_U);

	void NonDimToDimcalc(BoxData<double,NUMCOMPS>& a_U);

	void step_spherical_2O(LevelBoxData<double,NUMCOMPS>& a_Rhs,
			  LevelBoxData<double,NUMCOMPS>& a_U,
			  MHDLevelDataState& a_State,
			  double& a_min_dt);

	
	void step(LevelBoxData<double,NUMCOMPS>& a_Rhs,
			  LevelBoxData<double,NUMCOMPS>& a_JU_ave,
			  MHDLevelDataState& a_State,
			  double& a_min_dt);

	void step_spherical(LevelBoxData<double,NUMCOMPS>& a_Rhs,
			  LevelBoxData<double,NUMCOMPS>& a_JU_ave,
			  MHDLevelDataState& a_State,
			  double& a_min_dt);

	void consToPrimcalc(BoxData<double,NUMCOMPS>& a_W,
	                    const BoxData<double,NUMCOMPS>& a_U,
	                    const double gamma);

	void Fix_negative_P(BoxData<double,NUMCOMPS>& a_U,
	                    const double gamma);					

	void consToPrimSphcalc(BoxData<double,NUMCOMPS>& W_bar,
	                    const BoxData<double,NUMCOMPS>& a_U_sph,
	                    const BoxData<double,NUMCOMPS>& a_U_sph_actual,
	                    const double gamma);

	// If one is so motivated, one can expose internal pointwise functions.
	// It is not necessary for correctness, these are all internal functions
	// used by the "step" function above.
	PROTO_KERNEL_START
	void consToPrimF(State&         a_W,
	                 const State&   a_U,
	                 double a_gamma);

}
