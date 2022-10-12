#pragma once

#include "Proto.H"

//TODO: These variables need to go into a common file or be defined in a more accessible way
constexpr unsigned int NUMCOMPS=2*DIM+2;

#define NGHOST 6

using namespace Proto;

typedef Var<double,DIM> V;
typedef Var<double,NUMCOMPS> State;
class MHDLevelDataState;
namespace MHDOp
{


	// Compute the right-hand side for the MHD equations (i.e. div(flux)).
	// Note: a_Rhs is NOT multiplied by the time step, which may be required
	// by the RK4 implementation.

	// void step_spherical_2O(BoxData<double,NUMCOMPS>& a_Rhs,
	// 		               BoxData<double,NUMCOMPS>& a_Rhs_divB,
	// 					   double& a_min_dt,	
	//           const BoxData<double,NUMCOMPS>& a_U,
	//           const Box& a_rangeBox,
	// 		  BoxData<double,DIM>& a_x_sph_cc,
	// 		  BoxData<double,DIM>& a_x_sph_fc_1,
	// 		  BoxData<double,DIM>& a_x_sph_fc_2,
	// 		  BoxData<double,DIM>& a_x_sph_fc_3,
	// 		  BoxData<double,DIM>& a_dx_sph,
	// 		  BoxData<double,DIM>& a_face_area,
	// 		  BoxData<double,1>& a_cell_volume,
	//           const double a_dx,
	//           const double a_dy,
	//           const double a_dz,
	//           const double a_gamma,
	//           bool a_divB_calculated,
	//           bool a_min_dt_calculated);

	void step_spherical_2O(LevelBoxData<double,NUMCOMPS>& a_Rhs,
			  LevelBoxData<double,NUMCOMPS>& a_JU_ave,
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

	void consToPrimcalc(BoxData<double,NUMCOMPS>& W_bar,
	                    const BoxData<double,NUMCOMPS>& a_U_demapped,
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

	PROTO_KERNEL_START
	void upwindStateF(State& a_out,
	                  const State& a_low,
	                  const State& a_high,
	                  int a_dir,
	                  double a_gamma);


}
