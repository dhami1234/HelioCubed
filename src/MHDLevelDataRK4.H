#ifndef _MHD_LEVEL_DATA_RK4_H_
#define _MHD_LEVEL_DATA_RK4_H_

#include "Proto.H"
#include "Proto_LevelBoxData.H"
#include "MHDOp.H"
#include "MHD_Artificial_Viscosity.H"
#include "MHD_divB_Cleaning.H"
#include <iostream>
using namespace Proto;

#define NCOMPS 2*DIM+2

class MHDLevelDataDX;

class MHDLevelDataState
{
public:
	MHDLevelDataState();
	~MHDLevelDataState();
	MHDLevelDataState(const ProblemDomain& a_probDom,
	                  const Point& a_boxSize,
	                  const double a_dx,
	                  const double a_dy,
	                  const double a_dz,
	                  const double a_gamma);
	double m_dx;
	double m_dy;
	double m_dz;
	double m_gamma;
	bool m_divB_calculated;
	bool m_Viscosity_calculated;
	bool m_min_dt_calculated;
	double m_min_dt;
//TODO: Put in MEMTYPE so we can run on the GPU?
	LevelBoxData<double,NUMCOMPS> m_U;
	LevelBoxData<double,NUMCOMPS> m_U_old;
	LevelBoxData<double,NUMCOMPS> m_divB;
	LevelBoxData<double,NUMCOMPS> m_Viscosity;
	LevelBoxData<double,NUMCOMPS> m_BC;
	LevelBoxData<double,1> m_Jacobian_ave;
	LevelBoxData<double,DIM*DIM> m_N_ave_f;
	std::array<LevelBoxData<double,DIM>,DIM> m_NT;
	LevelBoxData<double,DIM> m_X_corners;
	LevelBoxData<double> m_J;

// Things needed by Spherical mapping with radial flow correction
	LevelBoxData<double,DIM*DIM> m_A_1_avg;
	LevelBoxData<double,DIM*DIM> m_A_2_avg;
	LevelBoxData<double,DIM*DIM> m_A_3_avg;
	LevelBoxData<double,DIM*DIM> m_A_inv_1_avg;
	LevelBoxData<double,DIM*DIM> m_A_inv_2_avg;
	LevelBoxData<double,DIM*DIM> m_A_inv_3_avg;
	LevelBoxData<double,DIM*DIM> m_detAA_avg;
	LevelBoxData<double,DIM*DIM> m_detAA_inv_avg;
	LevelBoxData<double,1> m_r2rdot_avg;
	LevelBoxData<double,1> m_detA_avg;
	LevelBoxData<double,DIM> m_A_row_mag_avg;
	LevelBoxData<double,1> m_r2detA_1_avg;
	LevelBoxData<double,DIM*DIM> m_r2detAA_1_avg;
	LevelBoxData<double,DIM> m_r2detAn_1_avg;
	LevelBoxData<double,DIM> m_A_row_mag_1_avg;
	LevelBoxData<double,1> m_rrdotdetA_2_avg;
	LevelBoxData<double,DIM*DIM> m_rrdotdetAA_2_avg;
	LevelBoxData<double,DIM> m_rrdotd3ncn_2_avg;
	LevelBoxData<double,DIM> m_A_row_mag_2_avg;
	LevelBoxData<double,1> m_rrdotdetA_3_avg;
	LevelBoxData<double,DIM*DIM> m_rrdotdetAA_3_avg;
	LevelBoxData<double,DIM> m_rrdotncd2n_3_avg;
	LevelBoxData<double,DIM> m_A_row_mag_3_avg;
// Things needed by Second order Spherical
	LevelBoxData<double,DIM> m_x_sph_cc;
	LevelBoxData<double,DIM> m_x_sph_fc_1;
	LevelBoxData<double,DIM> m_x_sph_fc_2;
	LevelBoxData<double,DIM> m_x_sph_fc_3;
	LevelBoxData<double,DIM> m_dx_sph;
	LevelBoxData<double,DIM> m_face_area;
	LevelBoxData<double,1> m_cell_volume;

	DisjointBoxLayout m_dbl;
	ProblemDomain m_probDom;
	void increment(const MHDLevelDataDX& a_DX);
};

class MHDLevelDataDX
{
public:
	MHDLevelDataDX();
	~MHDLevelDataDX();
	void init(MHDLevelDataState& a_State);
	void increment(double& a_weight, const MHDLevelDataDX& a_incr);
	void operator*=(const double& a_weight);
//TODO: Put in MEMTYPE so we can run on the GPU?
	LevelBoxData<double,NUMCOMPS> m_DU;
	DisjointBoxLayout m_dbl;
};

class MHDLevelDataRK4Op
{
public:
	MHDLevelDataRK4Op();
	~MHDLevelDataRK4Op();
	void operator()(MHDLevelDataDX& a_DX,
	                double a_time,
	                double a_dt,
	                MHDLevelDataState& a_State);
	PROTO_KERNEL_START
	void fill_ghostsF(State& a_U,
	                  const double a_gamma);
};

class MHDLevelDatadivBOp
{
public:
	MHDLevelDatadivBOp();
	~MHDLevelDatadivBOp();
	void operator()(MHDLevelDataDX& a_DX,
	                double a_time,
	                double a_dt,
	                MHDLevelDataState& a_State);
};

class MHDLevelDataViscosityOp
{
public:
	MHDLevelDataViscosityOp();
	~MHDLevelDataViscosityOp();
	void operator()(MHDLevelDataDX& a_DX,
	                double a_time,
	                double a_dt,
	                MHDLevelDataState& a_State);
};

#endif
