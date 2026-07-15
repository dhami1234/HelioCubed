#pragma once

#include "Proto.H"
#include "MHDOp.H"

typedef Var<double,NUMCOMPS> State;

namespace MHD_Riemann_Solvers {
	void Physical_Flux(BoxData<double,NUMCOMPS>& a_F,
	                   const BoxData<double,NUMCOMPS>& a_W,
	                   const int a_d,
	                   const double a_gamma);

	void Rusanov_Solver(BoxData<double,NUMCOMPS>& a_F_f,
	                    const BoxData<double,NUMCOMPS>& a_W_low,
	                    const BoxData<double,NUMCOMPS>& a_W_high,
	                    const int a_d,
	                    const double a_gamma);

	void Roe8Wave_Solver(BoxData<double,NUMCOMPS>& a_F_f,
	                     const BoxData<double,NUMCOMPS>& a_W_low,
	                     const BoxData<double,NUMCOMPS>& a_W_high,
	                     const int a_d,
	                     const double a_gamma);
}
