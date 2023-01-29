#pragma once

#include "Proto.H"
#include "MHDOp.H"

typedef Var<double,NUMCOMPS> State;

namespace MHD_Riemann_Solvers {
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

	void Spherical_Riemann_Solver(BoxData<double,NUMCOMPS>& a_F_ave_f,
	                              const BoxData<double,NUMCOMPS>& a_W_ave_low,
	                              const BoxData<double,NUMCOMPS>& a_W_ave_high,
								  const BoxData<double,NUMCOMPS>& a_W_ave_low_actual,
	                              const BoxData<double,NUMCOMPS>& a_W_ave_high_actual,
	                              const BoxData<double,1>& a_r2detA_1_avg,
	                              const BoxData<double,DIM*DIM>& a_r2detAA_1_avg,
	                              const BoxData<double,DIM>& a_r2detAn_1_avg,
	                              const BoxData<double,DIM>& a_n_1_avg,
	                              const BoxData<double,DIM*DIM>& a_A_1_avg,
	                              const BoxData<double,1>& a_rrdotdetA_2_avg,
	                              const BoxData<double,DIM*DIM>& a_rrdotdetAA_2_avg,
	                              const BoxData<double,DIM>& a_rrdotd3ncn_2_avg,
								  const BoxData<double,DIM*DIM>& a_A_2_avg,
	                              const BoxData<double,1>& a_rrdotdetA_3_avg,
	                              const BoxData<double,DIM*DIM>& a_rrdotdetAA_3_avg,
	                              const BoxData<double,DIM>& a_rrdotncd2n_3_avg,
								  const BoxData<double,DIM*DIM>& a_A_3_avg,
	                              const int a_d,
	                              const double a_gamma,
								  const double a_dx,
	                    		  const double a_dy,
	                    		  const double a_dz);

	void MHDSphericalFlux(
					BoxData<double,NUMCOMPS,MEM>& a_flux,                     
					const BoxData<double,NUMCOMPS,MEM>& a_prim4,
					const BoxData<double,NUMCOMPS,MEM>& a_prim2,
					const BoxData<double,NUMCOMPS,MEM>& a_prim_actual4,
					const BoxData<double,NUMCOMPS,MEM>& a_prim_actual2,
					const BoxData<double,DIM,MEM,DIM>& a_DrDetAA4,
					const BoxData<double,DIM,MEM,DIM>& a_DrDetAA2,
					const BoxData<double,DIM,MEM,DIM>& a_A4,
					const BoxData<double,DIM,MEM,DIM>& a_A2,
					const BoxData<double,1,MEM>& a_DrDetA4,                    
					const BoxData<double,1,MEM>& a_DrDetA2,
					const BoxData<double,DIM,MEM>& a_DrAdjA4,                    
					const BoxData<double,DIM,MEM>& a_DrAdjA2,
					const double& a_gamma,
					int a_dir);
}
