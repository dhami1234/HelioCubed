#pragma once

#include "Proto.H"
#include "MHDOp.H"

typedef Var<double, NUMCOMPS> State;

namespace MHD_Riemann_Solvers
{
	void Rusanov_Solver(BoxData<double, NUMCOMPS> &a_F_f,
						const BoxData<double, NUMCOMPS> &a_W_low,
						const BoxData<double, NUMCOMPS> &a_W_high,
						const int a_d,
						const double a_gamma);

	void Roe8Wave_Solver(BoxData<double, NUMCOMPS> &a_F_f,
						 const BoxData<double, NUMCOMPS> &a_W_low,
						 const BoxData<double, NUMCOMPS> &a_W_high,
						 const int a_d,
						 const double a_gamma);

	void Spherical_Riemann_Solver(BoxData<double, NUMCOMPS> &a_F_ave_f,
								  const BoxData<double, NUMCOMPS> &a_W_ave_low,
								  const BoxData<double, NUMCOMPS> &a_W_ave_high,
								  const BoxData<double, NUMCOMPS> &a_W_ave_low_actual,
								  const BoxData<double, NUMCOMPS> &a_W_ave_high_actual,
								  const BoxData<double, 1> &a_Dr_detA_avg,
								  const BoxData<double, DIM * DIM> &a_Dr_detA_A_avg,
								  const BoxData<double, DIM> &a_Dr_adjA_avg,
								  const BoxData<double, DIM> &a_A_row_mag_face_avg,
								  const int a_d,
								  const double a_gamma,
								  const double a_dx,
								  const double a_dy,
								  const double a_dz,
								  const int a_order);

	void MHDSphericalFlux(
		BoxData<double, NUMCOMPS, MEM> &a_flux,
		const BoxData<double, NUMCOMPS, MEM> &a_prim4,
		const BoxData<double, NUMCOMPS, MEM> &a_prim2,
		const BoxData<double, NUMCOMPS, MEM> &a_prim_actual4,
		const BoxData<double, NUMCOMPS, MEM> &a_prim_actual2,
		const BoxData<double, DIM, MEM, DIM> &a_DrDetAA4,
		const BoxData<double, DIM, MEM, DIM> &a_DrDetAA2,
		const BoxData<double, 1, MEM> &a_DrDetA4,
		const BoxData<double, 1, MEM> &a_DrDetA2,
		const BoxData<double, DIM, MEM> &a_DrAdjA4,
		const BoxData<double, DIM, MEM> &a_DrAdjA2,
		const double &a_gamma,
		int a_dir);

	void MHDSphericalFlux_2O(
		BoxData<double, NUMCOMPS, MEM> &a_flux,
		const BoxData<double, NUMCOMPS, MEM> &a_prim4,
		const BoxData<double, NUMCOMPS, MEM> &a_prim2,
		const BoxData<double, NUMCOMPS, MEM> &a_prim_actual4,
		const BoxData<double, NUMCOMPS, MEM> &a_prim_actual2,
		const BoxData<double, DIM, MEM, DIM> &a_DrDetAA4,
		const BoxData<double, DIM, MEM, DIM> &a_DrDetAA2,
		const BoxData<double, 1, MEM> &a_DrDetA4,
		const BoxData<double, 1, MEM> &a_DrDetA2,
		const BoxData<double, DIM, MEM> &a_DrAdjA4,
		const BoxData<double, DIM, MEM> &a_DrAdjA2,
		const double &a_gamma,
		int a_dir);
}
