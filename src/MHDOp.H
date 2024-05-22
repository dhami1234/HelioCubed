#pragma once

#include "Proto.H"
#include "MHD_Input_Parsing.H"

extern Parsefrominputs inputs;

//TODO: These variables need to go into a common file or be defined in a more accessible way
constexpr unsigned int NUMCOMPS=2*DIM+2 + TURB*3;

#define NGHOST 2
#define CRHO 0
#define CVELSTART 1
#define CBSTART 5
#define CPRES 4
#define CENG 4

#define iRHO 0
#define iVX 1
#define iMOMX 1
#define iVY 2
#define iMOMY 2
#define iVZ 3
#define iMOMZ 3
#define iE 4
#define iP 4
#define iBX 5
#define iBY 6
#define iBZ 7
#if TURB == 1
	#define iRHOZ2 8
	#define iZ2 8
	#define iRHOZ2SIGMA 9
	#define iSIGMA 9
	#define iRHOLAMBDA 10
	#define iLAMBDA 10
#endif


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
