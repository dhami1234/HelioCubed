#ifndef _MHD_LEVEL_DATA_RK4_H_
#define _MHD_LEVEL_DATA_RK4_H_

#include "Proto.H"
#include "Proto_LevelBoxData.H"
#include "MHDOp.H"
#include <iostream>
using namespace Proto;

constexpr MemType MEM = MEMTYPE_DEFAULT;

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
	bool m_divV_calculated;
	bool m_Viscosity_calculated;
	bool m_min_dt_calculated;
	bool m_CME_inserted;
	bool m_CME_checkpoint_written;
	double m_min_dt;
	double m_time;
//TODO: Put in MEMTYPE so we can run on the GPU?
	LevelBoxData<double,NUMCOMPS> m_U;
	LevelBoxData<double,NUMCOMPS> m_U_old;
	LevelBoxData<double,NUMCOMPS> m_divB;
	LevelBoxData<double,1> m_divV;
	LevelBoxData<double,NUMCOMPS> m_Viscosity;
	LevelBoxData<double,NUMCOMPS> m_BC;
	LevelBoxData<double,NUMCOMPS> m_CME;
	LevelBoxData<double,1> m_Jacobian_ave;
	LevelBoxData<double,DIM*DIM> m_N_ave_f;
	Array<LevelBoxData<double,DIM>,DIM> m_NT;
	LevelBoxData<double,DIM> m_X_corners;
	LevelBoxData<double> m_J;

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

#endif
