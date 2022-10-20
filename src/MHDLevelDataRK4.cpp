#include "MHDLevelDataRK4.H"
#include "MHD_Set_Boundary_Values.H"
#include "MHD_Mapping.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Constants.H"
#include "PolarExchangeCopier.H"

// For Chrono Timer (Talwinder)
#include <chrono>
#include <iostream>
#include <iomanip>

extern Parsefrominputs inputs;

MHDLevelDataState::MHDLevelDataState()
{
}

MHDLevelDataState::~MHDLevelDataState()
{
}

MHDLevelDataState::MHDLevelDataState(const ProblemDomain& a_probDom,
                                     const Point& a_boxSize,
                                     const double a_dx,
                                     const double a_dy,
                                     const double a_dz,
                                     const double a_gamma) :
    m_dx(a_dx),
    m_dy(a_dy),
    m_dz(a_dz),
	m_gamma(a_gamma),
	m_dbl(a_probDom, a_boxSize),
	m_probDom(a_probDom),
	m_min_dt(0.0)
{
	m_U.define(m_dbl,Point::Zero());
	m_U_old.define(m_dbl,Point::Zero());
	m_divB.define(m_dbl,Point::Ones(NGHOST));
	m_Viscosity.define(m_dbl,Point::Ones(NGHOST));
	m_BC.define(m_dbl,Point::Ones(NGHOST));
	m_Jacobian_ave.define(m_dbl,Point::Ones(NGHOST));
	m_N_ave_f.define(m_dbl,Point::Ones(NGHOST));
	for (int i=0; i<DIM; i++){
		m_NT[i].define(m_dbl,Point::Ones(NGHOST));
	}
	m_X_corners.define(m_dbl,Point::Ones(NGHOST));
	m_J.define(m_dbl,Point::Ones(NGHOST));

    m_A_1_avg.define(m_dbl,Point::Ones(NGHOST));
    m_A_2_avg.define(m_dbl,Point::Ones(NGHOST));
    m_A_3_avg.define(m_dbl,Point::Ones(NGHOST));
	m_A_inv_1_avg.define(m_dbl,Point::Ones(NGHOST));
    m_A_inv_2_avg.define(m_dbl,Point::Ones(NGHOST));
    m_A_inv_3_avg.define(m_dbl,Point::Ones(NGHOST));
    m_detAA_avg.define(m_dbl,Point::Ones(NGHOST));
    m_detAA_inv_avg.define(m_dbl,Point::Ones(NGHOST));
    m_r2rdot_avg.define(m_dbl,Point::Ones(NGHOST));
    m_detA_avg.define(m_dbl,Point::Ones(NGHOST));
    m_A_row_mag_avg.define(m_dbl,Point::Ones(NGHOST));
    m_r2detA_1_avg.define(m_dbl,Point::Ones(NGHOST));
    m_r2detAA_1_avg.define(m_dbl,Point::Ones(NGHOST));
    m_r2detAn_1_avg.define(m_dbl,Point::Ones(NGHOST));
	m_A_row_mag_1_avg.define(m_dbl,Point::Ones(NGHOST));
    m_rrdotdetA_2_avg.define(m_dbl,Point::Ones(NGHOST));
    m_rrdotdetAA_2_avg.define(m_dbl,Point::Ones(NGHOST));
    m_rrdotd3ncn_2_avg.define(m_dbl,Point::Ones(NGHOST));
	m_A_row_mag_2_avg.define(m_dbl,Point::Ones(NGHOST));
    m_rrdotdetA_3_avg.define(m_dbl,Point::Ones(NGHOST));
    m_rrdotdetAA_3_avg.define(m_dbl,Point::Ones(NGHOST));
    m_rrdotncd2n_3_avg.define(m_dbl,Point::Ones(NGHOST));
	m_A_row_mag_3_avg.define(m_dbl,Point::Ones(NGHOST));

	m_x_sph_cc.define(m_dbl,Point::Ones(NGHOST));
	m_x_sph_fc_1.define(m_dbl,Point::Ones(NGHOST));
	m_x_sph_fc_2.define(m_dbl,Point::Ones(NGHOST));
	m_x_sph_fc_3.define(m_dbl,Point::Ones(NGHOST));
	m_dx_sph.define(m_dbl,Point::Ones(NGHOST));
	m_face_area.define(m_dbl,Point::Ones(NGHOST));
	m_cell_volume.define(m_dbl,Point::Ones(NGHOST));
}

void MHDLevelDataState::increment(const MHDLevelDataDX& a_DX)
{
	m_U_old.setToZero();
	for (auto dit : a_DX.m_DU){
		m_U_old[ dit]+=m_U[ dit]; // This is needed in artificial viscosity calculation
		m_U[ dit]+=a_DX.m_DU[ dit];
	}
}

MHDLevelDataDX::MHDLevelDataDX()
{
}

MHDLevelDataDX::~MHDLevelDataDX()
{
}

void MHDLevelDataDX::init(MHDLevelDataState& a_State)
{
	m_dbl=a_State.m_dbl;
	m_DU.define(m_dbl,Point::Zero());
	// This leads to increasing iteration time for some reason. Treatment below does not. 
	// m_DU.setToZero();  // Update: Looks like it has been fixed now. // Maybe not
    for (auto dit : m_DU){
        m_DU[ dit].setVal(0.0);
    }
}

void MHDLevelDataDX::increment(double& a_weight, const MHDLevelDataDX& a_incr)
{
	for (auto dit : m_DU){	
		const BoxData<double,NUMCOMPS>& incr=a_incr.m_DU[ dit];
		BoxData<double,NUMCOMPS> temp(incr.box());
		incr.copyTo(temp);
		temp*=a_weight;
		m_DU[ dit]+=temp;
	}
}

void MHDLevelDataDX::operator*=(const double& a_weight)
{
	for (auto dit : m_DU){	
		m_DU[ dit]*=a_weight;
	}
}

MHDLevelDataRK4Op::MHDLevelDataRK4Op()
{
}

MHDLevelDataRK4Op::~MHDLevelDataRK4Op()
{
}

void MHDLevelDataRK4Op::operator()(MHDLevelDataDX& a_DX,
                                   double a_time,
                                   double a_dt,
                                   MHDLevelDataState& a_State
								   )
{	
	PR_TIME("MHDLevelDataRK4Op::operator");
	LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
	auto idOp = (1.0)*Shift(Point::Zeros());
	// (a_State.m_U).copyTo(new_state); // LevelBoxData copyTo doesn't copy ghost cells. Needs exchange if this is used
	for (auto dit : new_state){	
		(a_State.m_U[ dit]).copyTo(new_state[ dit]);
		// new_state[ dit]+=(a_DX.m_DU)[ dit];
        new_state[ dit]+=idOp((a_DX.m_DU)[ dit]);  //Phil found doing this improves performance
	}
	if (inputs.grid_type_global == 2){
		new_state.defineExchange<PolarExchangeCopier>(2,1);
	}
	new_state.exchange(); 

	
	
	if (inputs.LowBoundType != 0 || inputs.HighBoundType != 0) {
		if (inputs.Spherical_2nd_order == 0){
			MHD_Set_Boundary_Values::Set_Boundary_Values(new_state,a_State);			
		}
		if (inputs.Spherical_2nd_order == 1){
			MHD_Set_Boundary_Values::Set_Boundary_Values_Spherical_2O(new_state,a_State);
		}
	} 
	
	double dt_temp = 1.0e10;
	double dt_new;
		
	
	if (inputs.grid_type_global == 2){
		if (inputs.Spherical_2nd_order == 0){
			MHDOp::step_spherical(a_DX.m_DU,new_state,a_State, dt_temp);
		}
		if (inputs.Spherical_2nd_order == 1){
			PR_TIME("MHDOp::step_spherical_2O");
			MHDOp::step_spherical_2O(a_DX.m_DU,new_state,a_State, dt_temp);
		}
	} else {
		MHDOp::step(a_DX.m_DU,new_state,a_State, dt_temp);
	}

	if (!a_State.m_min_dt_calculated){
		double mintime;
		#ifdef PR_MPI
			MPI_Reduce(&dt_temp, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
			MPI_Bcast(&mintime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif
		a_State.m_min_dt = mintime;
	}
	a_State.m_divB_calculated = true; // This makes sure that divB is calculated only once in RK4
	a_State.m_Viscosity_calculated = true; // This makes sure that Viscosity is calculated only once in RK4
	a_State.m_min_dt_calculated = true; // This makes sure that min_dt is calculated only once in RK4
	a_DX*=a_dt;
}



MHDLevelDatadivBOp::MHDLevelDatadivBOp()
{
}

MHDLevelDatadivBOp::~MHDLevelDatadivBOp()
{
}

void MHDLevelDatadivBOp::operator()(MHDLevelDataDX& a_DX,
                                    double a_time,
                                    double a_dt,
                                    MHDLevelDataState& a_State)
{
	a_State.m_divB.copyTo(a_DX.m_DU);
	a_DX*=a_dt;
}

MHDLevelDataViscosityOp::MHDLevelDataViscosityOp()
{
}

MHDLevelDataViscosityOp::~MHDLevelDataViscosityOp()
{
}

void MHDLevelDataViscosityOp::operator()(MHDLevelDataDX& a_DX,
                                         double a_time,
                                         double a_dt,
                                         MHDLevelDataState& a_State)
{
	LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
	auto idOp = (1.0)*Shift(Point::Zeros());
	// (a_State.m_U).copyTo(new_state); // LevelBoxData copyTo doesn't copy ghost cells. Needs exchange if this is used
	for (auto dit : new_state){	
		(a_State.m_U[ dit]).copyTo(new_state[ dit]);
		// (a_State.m_U_old[ dit]).copyTo(new_state[ dit]);
		// new_state[ dit]+=(a_DX.m_DU)[ dit];
        new_state[ dit]+=idOp((a_DX.m_DU)[ dit]);  //Phil found doing this improves performance
	}
	if (inputs.grid_type_global == 2){
		new_state.defineExchange<PolarExchangeCopier>(2,1);
	}
	new_state.exchange(); 

	
   if (inputs.LowBoundType != 0 || inputs.HighBoundType != 0) {
		if (inputs.Spherical_2nd_order == 0){
			MHD_Set_Boundary_Values::Set_Boundary_Values(new_state,a_State);			
		}
		if (inputs.Spherical_2nd_order == 1){
			MHD_Set_Boundary_Values::Set_Boundary_Values_Spherical_2O(new_state,a_State);
		}
	}

	if (inputs.grid_type_global == 2){
		if (inputs.Spherical_2nd_order == 0){
			MHD_Artificial_Viscosity::step_spherical(a_DX.m_DU,new_state,a_State);
		}
	} else {
		MHD_Artificial_Viscosity::step(a_DX.m_DU,new_state,a_State);
	}

	a_DX*=a_dt;
}


