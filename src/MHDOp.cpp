#include "Proto.H"
#include "MHDOp.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
// #include "Proto_WriteBoxData.H"
// For Chrono Timer (Talwinder)
#include <chrono>
#include <iostream>
#include <iomanip>

//////////////////////////////
#include "MHD_Limiters.H"
#include "MHD_Mapping.H"
#include "MHD_Riemann_Solvers.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Constants.H"
#include "MHD_CFL.H"
#include "MHDLevelDataRK4.H"

extern Parsefrominputs inputs;

typedef BoxData<double,1,HOST> Scalar;
typedef BoxData<double,NUMCOMPS,HOST> Vector;
/// @brief MHDOp namespace
namespace MHDOp {
	/**
	 * @brief Function to covert conserevd variables to primitive variables. 
	 * @param a_W the output primitive variables.
	 * @param a_U the input conserved variables.
	 * @param a_gamma gamma.
	 */ 
	PROTO_KERNEL_START
	void
	consToPrimF(State&         a_W,
	            const State&   a_U,
	            double a_gamma)
	{
		double rho = a_U(0);
		double v2 = 0.0;
		double B2 = 0.0;
		double gamma = a_gamma;
		a_W(0) = rho;

		for (int i = 1; i <= DIM; i++)
		{
			double v, B;
			v = a_U(i) / rho;
			B = a_U(DIM+1+i);
			a_W(i) = v;
			a_W(DIM+1+i) = a_U(DIM+1+i);
			v2 += v*v;
			B2 += B*B;
		}

		a_W(NUMCOMPS-1-DIM) = (a_U(NUMCOMPS-1-DIM) - .5 * rho * v2  - B2/8.0/c_PI) * (gamma - 1.0);

	}
	PROTO_KERNEL_END(consToPrimF, consToPrim)


	/**
	 * @brief Function to covert conserevd variables to primitive variables. 
	 * @param a_W the output primitive variables.
	 * @param a_U the input conserved variables.
	 * @param a_gamma gamma.
	 */ 
	void consToPrimcalc(BoxData<double,NUMCOMPS>& a_W,
	                    const BoxData<double,NUMCOMPS>& a_U,
	                    const double gamma)
	{
		a_W = forall<double,NUMCOMPS>(consToPrim,a_U, gamma);
	}


	/**
	 * @brief Function to covert dimensional variables to non-dimensional variables. 
	 * @param a_U will get converted from dimensional to non-dimensional.
	 */  
	PROTO_KERNEL_START
	void
	DimToNonDimF(const Point& a_pt,
				 State& a_U)
	{
		// When going from cgs to non-dimensional
		// density is divided by a scale density density_scale
		// velocity is divided by a scale velocity velocity_scale
		// energy is divided by density_scale*velocity_scale*velocity_scale
		// magnetic field is divided by sqrt(density_scale*velocity_scale*velocity_scale)
		// pref = (lismN*eos_mp)*lismV*lismV	
		// CMEData(CHF_IX[i;j;k],0) = (mCME)/(lismN*eos_mp) !rhoTD
		// CMEData(CHF_IX[i;j;k],1) = B_x_dom/sqrt( pref ) !Bx
		// CMEData(CHF_IX[i;j;k],2) = B_y_dom/sqrt( pref ) !By
		// CMEData(CHF_IX[i;j;k],3) = B_z_dom/sqrt( pref ) !Bz
		// CMEData(CHF_IX[i;j;k],4) = energy_control*e0TD/pref !e0TD
		// CMEData(CHF_IX[i;j;k],5) = V_x_dom*1.0e5/lismV !Vx
		// CMEData(CHF_IX[i;j;k],6) = V_y_dom*1.0e5/lismV !Vy
		// CMEData(CHF_IX[i;j;k],7) = V_z_dom*1.0e5/lismV !Vz
		double pref = (inputs.density_scale*c_MP*inputs.velocity_scale*inputs.velocity_scale);
		#if DIM == 2
		a_U(0) = a_U(0)/(inputs.density_scale*c_MP);
		a_U(1) = a_U(1)/(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(2) = a_U(2)/(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(3) = a_U(3)/pref;
		a_U(4) = a_U(4)/sqrt(pref);
		a_U(5) = a_U(5)/sqrt(pref);
		#endif

		#if DIM == 3
		a_U(0) = a_U(0)/(inputs.density_scale*c_MP);
		a_U(1) = a_U(1)/(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(2) = a_U(2)/(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(3) = a_U(3)/(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(4) = a_U(4)/pref;
		a_U(5) = a_U(5)/sqrt(pref);
		a_U(6) = a_U(6)/sqrt(pref);
		a_U(7) = a_U(7)/sqrt(pref);
		#endif
	}
	PROTO_KERNEL_END(DimToNonDimF, DimToNonDim)


	/**
	 * @brief Function to covert dimensional variables to non-dimensional variables. 
	 * @param a_U will get converted from dimensional to non-dimensional.
	 */ 
	void DimToNonDimcalc(BoxData<double,NUMCOMPS>& a_U)
	{
		forallInPlace_p(DimToNonDim, a_U);
	}

	/**
	 * @brief Function to covert non-dimensional variables to dimensional variables. 
	 * @param a_U will get converted from non-dimensional to dimensional.
	 */  
	PROTO_KERNEL_START
	void
	NonDimToDimF(const Point& a_pt,
				 State& a_U)
	{
		// When going from non-dimensional to cgs
		// density is multiplied by a scale density density_scale
		// velocity is multiplied by a scale velocity velocity_scale
		// energy is multiplied by density_scale*velocity_scale*velocity_scale
		// magnetic field is multiplied by sqrt(density_scale*velocity_scale*velocity_scale)

		double pref = (inputs.density_scale*c_MP*inputs.velocity_scale*inputs.velocity_scale);
		#if DIM == 2
		a_U(0) = a_U(0)*(inputs.density_scale*c_MP);
		a_U(1) = a_U(1)*(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(2) = a_U(2)*(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(3) = a_U(3)*pref;
		a_U(4) = a_U(4)*sqrt(pref);
		a_U(5) = a_U(5)*sqrt(pref);
		#endif

		#if DIM == 3
		a_U(0) = a_U(0)*(inputs.density_scale*c_MP);
		a_U(1) = a_U(1)*(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(2) = a_U(2)*(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(3) = a_U(3)*(inputs.density_scale*c_MP*inputs.velocity_scale);
		a_U(4) = a_U(4)*pref;
		a_U(5) = a_U(5)*sqrt(pref);
		a_U(6) = a_U(6)*sqrt(pref);
		a_U(7) = a_U(7)*sqrt(pref);
		#endif
	}
	PROTO_KERNEL_END(NonDimToDimF, NonDimToDim)


	/**
	 * @brief Function to covert non-dimensional variables to dimensional variables. 
	 * @param a_U will get converted from non-dimensional to dimensional.
	 */ 
	void NonDimToDimcalc(BoxData<double,NUMCOMPS>& a_U)
	{
		forallInPlace_p(NonDimToDim, a_U);
	}


	/**
	 * @brief Function to covert conserevd variables to primitive variables, used in spherical mapping. 
	 * @param a_W_sph the output primitive variables.
	 * @param a_U_sph the input conserved variables, but the v and B vectors are scaled with row magnitudes of A matrix (see overleaf document).
	 * @param a_U_sph_actual the input conserved variables.
	 * @param a_gamma gamma.
	 */
	PROTO_KERNEL_START
	void
	consToPrimSphF(State&         a_W_sph,
	            const State&   a_U_sph,
	            const State&   a_U_sph_actual,
	            double a_gamma)
	{
		double rho = a_U_sph(0);
		double v2 = 0.0;
		double B2 = 0.0;
		double gamma = a_gamma;
		a_W_sph(0) = rho;

		for (int i = 1; i <= DIM; i++)
		{
			double v, v_actual, B, B_actual;
			v = a_U_sph(i) / rho;
			v_actual = a_U_sph_actual(i) / rho;
			B = a_U_sph(DIM+1+i);
			B_actual = a_U_sph_actual(DIM+1+i);
			a_W_sph(i) = v;
			a_W_sph(DIM+1+i) = a_U_sph(DIM+1+i);
			v2 += v_actual*v_actual;
			B2 += B_actual*B_actual;
		}

		a_W_sph(NUMCOMPS-1-DIM) = (a_U_sph(NUMCOMPS-1-DIM) - .5 * rho * v2  - B2/8.0/c_PI) * (gamma - 1.0);
		// a_W_sph(NUMCOMPS-1-DIM) = a_U_sph(NUMCOMPS-1-DIM);

	}
	PROTO_KERNEL_END(consToPrimSphF, consToPrimSph)

	/**
	 * @brief Function to covert conserevd variables to primitive variables, used in spherical mapping. 
	 * @param a_W_sph the output primitive variables.
	 * @param a_U_sph the input conserved variables, but the v and B vectors are scaled with row magnitudes of A matrix (see overleaf document).
	 * @param a_U_sph_actual the input conserved variables.
	 * @param a_gamma gamma.
	 */
	void consToPrimSphcalc(BoxData<double,NUMCOMPS>& a_W_sph,
	                    const BoxData<double,NUMCOMPS>& a_U_sph,
	                    const BoxData<double,NUMCOMPS>& a_U_sph_actual,
	                    const double gamma)
	{
		a_W_sph = forall<double,NUMCOMPS>(consToPrimSph,a_U_sph,a_U_sph_actual, gamma);
	}

	/**
	 * @brief Function to calculate the Powell term, without the divB multiplied to it. 
	 * @param a_P the output terms.
	 * @param a_W the input primitive variables.
	 */
	PROTO_KERNEL_START
	void PowellF(State&         a_P,
	             const State&   a_W)
	{

#if DIM==2
		a_P(0) = 0.;
		a_P(1) = a_W(4)/4.0/c_PI;
		a_P(2) = a_W(5)/4.0/c_PI;
		a_P(3) = a_W(1)*a_W(4)/4.0/c_PI + a_W(2)*a_W(5)/4.0/c_PI;
		a_P(4) = a_W(1);
		a_P(5) = a_W(2);
#endif

#if DIM==3
		a_P(0) = 0.;
		a_P(1) = a_W(5)/4.0/c_PI;
		a_P(2) = a_W(6)/4.0/c_PI;
		a_P(3) = a_W(7)/4.0/c_PI;
		a_P(4) = a_W(1)*a_W(5)/4.0/c_PI + a_W(2)*a_W(6)/4.0/c_PI + a_W(3)*a_W(7)/4.0/c_PI;
		a_P(5) = a_W(1);
		a_P(6) = a_W(2);
		a_P(7) = a_W(3);

#endif
	}
	PROTO_KERNEL_END(PowellF, Powell)


	/**
	 * @brief Function to average cell centered B on both sides of a face. 
	 * @param a_Bavg the output BoxData.
	 * @param a_W_ave the input cell averaged BoxData of primitive variables.
	 * @param a_d direction to average.
	 */ 
	PROTO_KERNEL_START
	void BavgcalcF(State& a_Bavg,
	               const State& a_W_ave,
	               int a_d)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			a_Bavg(i) = a_W_ave(2+DIM+a_d);
		}
	}
	PROTO_KERNEL_END(BavgcalcF, Bavgcalc)
	
	/**
	 * @brief Function to transfer data from lower dimensional BoxData to a higher dimensional BoxData. 
	 * @param a_F the output BoxData.
	 * @param a_F_temp the input BoxData.
	 * @param a_s where to transfer.
	 */ 
	PROTO_KERNEL_START
	void Fill_flux_calcF(const Point& a_pt,
						Var<double,DIM, MEM ,NUMCOMPS>&       a_F,
						const Var<double,NUMCOMPS>&        a_F_temp,
	            		const int a_s)
	{
		for (int i=0; i<NUMCOMPS; i++){
			a_F(a_s,i) = a_F_temp(i);
		}
	}
	PROTO_KERNEL_END(Fill_flux_calcF, Fill_flux_calc)

	/**
	 * @brief Function to transpose BoxData. 
	 * @param a_F the output BoxData.
	 * @param a_F_temp the input BoxData.
	 */ 
	PROTO_KERNEL_START
	void Transpose_calcF(const Point& a_pt,
						Var<double,NUMCOMPS>&       a_F,
						const Var<double,1,MEM,NUMCOMPS>&        a_F_temp)
	{
		for (int i=0; i<NUMCOMPS; i++){
			a_F(i) = a_F_temp(0,i);
		}
	}
	PROTO_KERNEL_END(Transpose_calcF, Transpose_calc)

	
	/**
	 * @brief Function to calculate the right hand side of the finite volume method, Magnetic field divergence, and min dt according to CFL condition. 
	 * @param a_Rhs the output LevelBoxData.
	 * @param a_JU_ave the input LevelBoxData containing 4th order averaged product of Jacobian and conserved variables.
	 * @note If no mapping is used, J = 1.
	 */ 
	void step(LevelBoxData<double,NUMCOMPS>& a_Rhs,
			  LevelBoxData<double,NUMCOMPS>& a_JU_ave,
			  MHDLevelDataState& a_State,
			  double& a_min_dt)
	{

		
		static Stencil<double> m_laplacian;
		static Stencil<double> m_deconvolve;
		static Stencil<double> m_copy;
		static Stencil<double> m_laplacian_f[DIM];
		static Stencil<double> m_deconvolve_f[DIM];
		static Stencil<double> m_convolve_f[DIM];
		static Stencil<double> m_interp_H[DIM];
		static Stencil<double> m_interp_L[DIM];
		static Stencil<double> m_interp_edge[DIM];
		static Stencil<double> m_divergence[DIM];
		static Stencil<double> m_derivative[DIM];
		static bool initialized = false;
		if(!initialized)
		{
			m_laplacian = Stencil<double>::Laplacian();
			m_deconvolve = (-1.0/24.0)*m_laplacian + (1.0)*Shift(Point::Zeros());
			m_copy = 1.0*Shift(Point::Zeros());
			for (int dir = 0; dir < DIM; dir++)
			{
				m_laplacian_f[dir] = Stencil<double>::LaplacianFace(dir);
				m_deconvolve_f[dir] = (-1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
				m_convolve_f[dir] = (1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
				m_interp_H[dir] = Stencil<double>::CellToFaceH(dir);
				m_interp_L[dir] = Stencil<double>::CellToFaceL(dir);
				m_interp_edge[dir] = Stencil<double>::CellToFace(dir);
				m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
				m_derivative[dir] = Stencil<double>::Derivative(1,dir,2);
			}
			initialized =  true;
		}



		using namespace std;
		double a_dx = a_State.m_dx;
		double a_dy = a_State.m_dy;
		double a_dz = a_State.m_dz;
		if (DIM == 2) a_dz = 1.0;
		double volume = a_dx*a_dy*a_dz;
		double gamma = a_State.m_gamma;
		double dxd[3] = {a_dx, a_dy, a_dz};
		double dt_new;
		for (auto dit : a_State.m_U){
			Box dbx0 = a_JU_ave[dit].box();
			Box dbx1 = dbx0;
			Box dbx2 = dbx0.grow(0-NGHOST);

			a_Rhs[dit].setVal(0.0);
			if (!a_State.m_divB_calculated) a_State.m_divB[dit].setVal(0.0);

			HDF5Handler h5;

			auto a_U_ave = Operator::_cellTensorQuotient(a_JU_ave[dit],a_State.m_J[dit],a_JU_ave[dit],a_State.m_J[dit]);

			Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U_ave, gamma);
			Vector U = m_deconvolve(a_U_ave);
			Vector W  = forall<double,NUMCOMPS>(consToPrim,U, gamma);
			Vector W_ave = m_laplacian(W_bar,1.0/24.0);
			W_ave += W;
			if (!a_State.m_min_dt_calculated){ 
				MHD_CFL::Min_dt_calc_func(dt_new, W_ave, dbx0, a_dx, a_dy, a_dz, gamma);	
				if (dt_new < a_min_dt) a_min_dt = dt_new;
			}

			for (int d = 0; d < DIM; d++)
			{
				Vector W_ave_low_temp(dbx0), W_ave_high_temp(dbx0);
				Vector W_ave_low(dbx0), W_ave_high(dbx0);
				W_ave_low_temp = m_interp_L[d](W_ave);
				W_ave_high_temp = m_interp_H[d](W_ave);
				Vector W_ave_edge = m_interp_edge[d](W_ave);
				MHD_Limiters::MHD_Limiters_4O(W_ave_low,W_ave_high,W_ave_low_temp,W_ave_high_temp,W_ave,W_bar,d,a_dx, a_dy, a_dz);
				Vector W_low = m_deconvolve_f[d](W_ave_low);
				Vector W_high = m_deconvolve_f[d](W_ave_high);
				BoxData<double,DIM,MEM,NUMCOMPS> F_f(dbx1);
				BoxData<double,DIM,MEM,NUMCOMPS> F_ave_f(dbx1);
				BoxData<double,NUMCOMPS> F_f_temp(dbx1);
				Vector F_f_mapped(dbx1);
				F_f_mapped.setVal(0.0);
				BoxData<double,DIM,MEM,NUMCOMPS> BF_f(dbx1), BF_ave_f(dbx1);
				Vector BF_f_mapped(dbx1);
				if (!a_State.m_divB_calculated) {BF_f_mapped.setVal(0.0);}
				double dx_d = dxd[d];
				for (int s = 0; s < DIM; s++) {
					if (inputs.Riemann_solver_type == 1) {
						MHD_Riemann_Solvers::Rusanov_Solver(F_f_temp,W_low,W_high,s,gamma);
					}
					if (inputs.Riemann_solver_type == 2) {
						MHD_Riemann_Solvers::Roe8Wave_Solver(F_f_temp,W_low,W_high,s,gamma);
					}
					forallInPlace_p(Fill_flux_calc, F_f, F_f_temp, s);
					if (!a_State.m_divB_calculated){
						Vector B_ave_f = forall<double,NUMCOMPS>(Bavgcalc, W_ave_edge, s);
						forallInPlace_p(Fill_flux_calc, BF_ave_f, B_ave_f, s);
					}
				}

				F_ave_f = m_convolve_f[d](F_f);
				BoxData<double,1,MEM,NUMCOMPS> fluxdir = Operator::_faceMatrixProductATB(a_State.m_NT[d][dit],F_ave_f,a_State.m_NT[d][dit],F_ave_f,d);
				forallInPlace_p(Transpose_calc, F_f_mapped, fluxdir);
				Vector Rhs_d = m_divergence[d](F_f_mapped);
				// Rhs_d *= -1./pow(dx_d,DIM);
				Rhs_d *= -1./volume;
				a_Rhs[dit] += Rhs_d;

				if (!a_State.m_divB_calculated){
					BoxData<double,1,MEM,NUMCOMPS> Bfluxdir = Operator::_faceMatrixProductATB(a_State.m_NT[d][dit],BF_ave_f,a_State.m_NT[d][dit],BF_ave_f,d);
					forallInPlace_p(Transpose_calc, BF_f_mapped, Bfluxdir);
					Vector B_Rhs_d = m_divergence[d](BF_f_mapped);
					// B_Rhs_d *= -1./pow(dx_d,DIM);
					B_Rhs_d *= -1./volume;
					a_State.m_divB[dit] += B_Rhs_d;
				}
			}
			if (!a_State.m_divB_calculated){
				Vector Powell_term = forall<double,NUMCOMPS>(Powell,W_ave);
				a_State.m_divB[dit] *= Powell_term;
			}
		}
	}



	// /**
	//  * @brief Function to calculate the right hand side of the finite volume method, Magnetic field divergence, and min dt according to CFL condition. 
	//  * This function is tailored for the special mapping of the spherical grid that preserves radial flow.
	//  * @param a_Rhs the output LevelBoxData.
	//  * @param a_JU_ave the input LevelBoxData containing 4th order averaged product of Jacobian and conserved variables.
	//  * @note If no mapping is used, J = 1.
	//  */ 
	// void step_spherical(LevelBoxData<double,NUMCOMPS>& a_Rhs,
	// 		  LevelBoxData<double,NUMCOMPS>& a_JU_ave,
	// 		  MHDLevelDataState& a_State,
	// 		  double& a_min_dt)
	// {	
	// 	static Stencil<double> m_laplacian;
	// 	static Stencil<double> m_deconvolve;
	// 	static Stencil<double> m_copy;
	// 	static Stencil<double> m_laplacian_f[DIM];
	// 	static Stencil<double> m_deconvolve_f[DIM];
	// 	static Stencil<double> m_convolve_f[DIM];
	// 	static Stencil<double> m_interp_H[DIM];
	// 	static Stencil<double> m_interp_L[DIM];
	// 	static Stencil<double> m_interp_edge[DIM];
	// 	static Stencil<double> m_interp_f_2nd[DIM];
	// 	static Stencil<double> m_divergence[DIM];
	// 	static bool initialized = false;
	// 	if(!initialized)
	// 	{
	// 		m_laplacian = Stencil<double>::Laplacian();
	// 		m_deconvolve = (-1.0/24.0)*m_laplacian + (1.0)*Shift(Point::Zeros());
	// 		m_copy = 1.0*Shift(Point::Zeros());
	// 		for (int dir = 0; dir < DIM; dir++)
	// 		{
	// 			m_laplacian_f[dir] = Stencil<double>::LaplacianFace(dir);
	// 			m_deconvolve_f[dir] = (-1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
	// 			m_convolve_f[dir] = (1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
	// 			m_interp_f_2nd[dir] = 0.5*Shift(Point::Zeros()) + 0.5*Shift(-Point::Basis(dir)); 
	// 			m_interp_H[dir] = Stencil<double>::CellToFaceH(dir);
	// 			m_interp_L[dir] = Stencil<double>::CellToFaceL(dir);
	// 			m_interp_edge[dir] = Stencil<double>::CellToFace(dir);
	// 			m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
	// 		}
	// 		initialized =  true;
	// 	}

	// 	using namespace std;
		
	// 	double a_dx = a_State.m_dx;
	// 	double a_dy = a_State.m_dy;
	// 	double a_dz = a_State.m_dz;
	// 	double gamma = a_State.m_gamma;
	// 	double dxd[3] = {a_dx, a_dy, a_dz}; // Because now its r, theta, phi
	// 	double dt_new;
	// 	for (auto dit : a_State.m_U){
	// 		Box dbx0 = a_JU_ave[dit].box();
	// 		Box dbx1 = dbx0.grow(NGHOST-NGHOST);
	// 		a_Rhs[dit].setVal(0.0);
	// 		Vector a_U_Sph_ave(dbx0), a_U_Sph_actual_ave(dbx0);
	// 		MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_ave, a_JU_ave[dit], a_State.m_detAA_inv_avg[dit], a_State.m_A_inv_avg[dit], a_State.m_r2rdot_avg[dit], a_State.m_detA_avg[dit], a_State.m_A_row_mag_avg[dit], false, 2);
	// 		MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_actual_ave, a_JU_ave[dit], a_State.m_detAA_inv_avg[dit], a_State.m_A_inv_avg[dit], a_State.m_r2rdot_avg[dit], a_State.m_detA_avg[dit], a_State.m_A_row_mag_avg[dit], true, 2);
	// 		MHD_Mapping::Correct_V_theta_phi_at_poles(a_U_Sph_ave, a_dx, a_dy, a_dz);
	// 		Vector W_bar = forall<double,NUMCOMPS>(consToPrimSph, a_U_Sph_ave, a_U_Sph_actual_ave, gamma);
	// 		Vector W_ave = m_copy(W_bar);
	// 		HDF5Handler h5;
	// 		// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, W_ave, "W_MHDOp");
	// 		// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, a_JU_ave[dit], "a_JU_ave_MHDOp");
	// 		if (!a_State.m_min_dt_calculated){ 
	// 			MHD_CFL::Min_dt_calc_func(dt_new, W_ave, dbx0, a_dx, a_dy, a_dz, gamma);	
	// 			if (dt_new < a_min_dt) a_min_dt = dt_new;
	// 		}
			
	// 		for (int d = 0; d < DIM; d++)
	// 		{
	// 			Vector W_ave_low_temp(dbx0), W_ave_high_temp(dbx0);
	// 			Vector W_ave_low(dbx0), W_ave_high(dbx0);
	// 			Vector W_ave_low_actual(dbx0), W_ave_high_actual(dbx0);

	// 			W_ave_low_temp = m_interp_L[d](W_ave);
	// 			W_ave_high_temp = m_interp_H[d](W_ave);

	// 			MHD_Limiters::MHD_Limiters_4O(W_ave_low,W_ave_high,W_ave_low_temp,W_ave_high_temp,W_ave,W_bar,d,a_dx, a_dy, a_dz);			
	// 			// MHD_Limiters::MHD_Limiters_minmod(W_ave_low,W_ave_high,W_ave,a_State.m_x_sph_cc[dit],a_State.m_dx_sph[dit],d);
	// 			MHD_Mapping::W_Sph_to_W_normalized_sph(W_ave_low_actual, W_ave_low, a_State.m_A_row_mag_1_avg[dit], a_State.m_A_row_mag_2_avg[dit], a_State.m_A_row_mag_3_avg[dit], d);
	// 			MHD_Mapping::W_Sph_to_W_normalized_sph(W_ave_high_actual, W_ave_low, a_State.m_A_row_mag_1_avg[dit], a_State.m_A_row_mag_2_avg[dit], a_State.m_A_row_mag_3_avg[dit], d);
	// 			// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, W_ave_low, "W_ave_low_bad" +to_string(d));
	// 			Vector F_ave_f(dbx0);
	// 			// Vector F_ave_f_sph(dbx0);
	// 			F_ave_f.setVal(0.0);
	// 			double dx_d = dxd[d];
	// 			MHD_Riemann_Solvers::Spherical_Riemann_Solver(F_ave_f, W_ave_low, W_ave_high, W_ave_low_actual, W_ave_high_actual, a_State.m_r2detA_1_avg[dit], a_State.m_r2detAA_1_avg[dit], a_State.m_r2detAn_1_avg[dit], a_State.m_n_1_avg[dit], a_State.m_A_1_avg[dit], a_State.m_rrdotdetA_2_avg[dit], a_State.m_rrdotdetAA_2_avg[dit], a_State.m_rrdotd3ncn_2_avg[dit],a_State.m_A_2_avg[dit], a_State.m_rrdotdetA_3_avg[dit], a_State.m_rrdotdetAA_3_avg[dit], a_State.m_rrdotncd2n_3_avg[dit],a_State.m_A_3_avg[dit], d, gamma, a_dx, a_dy, a_dz);	
	// 			// MHD_Mapping::JU_to_W_Sph_ave_calc_func(F_ave_f_sph, F_ave_f, (a_State.m_detAA_inv_avg)[ dit], (a_State.m_A_inv_avg)[ dit], (a_State.m_r2rdot_avg)[ dit], (a_State.m_detA_avg)[ dit], (a_State.m_A_row_mag_avg)[ dit], inputs.gamma, true);
	// 			if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, F_ave_f, "F_ave_f"+to_string(d));
	// 			// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, F_ave_f_sph, "F_ave_f_sph_bad"+to_string(d));
	// 			Vector Rhs_d = m_divergence[d](F_ave_f);
	// 			Rhs_d *= -1./dx_d;
	// 			a_Rhs[dit] += Rhs_d;
	// 		}
	// 		// Vector a_Rhs_sph(dbx0);
	// 		// MHD_Mapping::JU_to_W_Sph_ave_calc_func(a_Rhs_sph, a_Rhs[dit], (a_State.m_detAA_inv_avg)[ dit], (a_State.m_A_inv_avg)[ dit], (a_State.m_r2rdot_avg)[ dit], (a_State.m_detA_avg)[ dit], (a_State.m_A_row_mag_avg)[ dit], inputs.gamma, true);
	// 		// MHD_Mapping::JU_to_W_Sph_ave_calc_func2(a_Rhs_sph, a_Rhs[dit], (a_State.m_detAA_inv_avg)[ dit], (a_State.m_A_inv_avg)[ dit], (a_State.m_r2rdot_avg)[ dit], (a_State.m_detA_avg)[ dit], (a_State.m_A_row_mag_avg)[ dit], inputs.gamma, true);
	// 		// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, a_Rhs[dit], "a_Rhs_bad");
	// 		// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, a_Rhs_sph, "a_Rhs_sph_bad");
	// 	}
	// }



	template<class T,unsigned int CFLUX,unsigned int CPRIM, MemType MEM>
	BoxData<T,CFLUX,MEM>
	MHDSphericalFlux(                     
					const BoxData<T,CPRIM,MEM>& a_prim4,
					const BoxData<T,CPRIM,MEM>& a_prim2,
					const BoxData<T,CPRIM,MEM>& a_prim_actual4,
					const BoxData<T,CPRIM,MEM>& a_prim_actual2,
					const BoxData<T,DIM,MEM,DIM>& a_DrDetAA4,
					const BoxData<T,DIM,MEM,DIM>& a_DrDetAA2,
					const BoxData<T,DIM,MEM,DIM>& a_A4,
					const BoxData<T,DIM,MEM,DIM>& a_A2,
					const BoxData<T,1,MEM>& a_DrDetA4,                    
					const BoxData<T,1,MEM>& a_DrDetA2,
					const BoxData<T,DIM,MEM>& a_DrAdjA4,                    
					const BoxData<T,DIM,MEM>& a_DrAdjA2,
					const T& a_gamma,
					int a_dir)
	{
		#define CRHO 0
		#define CVELSTART 1
		#define CBSTART 5
		#define CPRES 4
		#define CENG 4
		
		auto wnorm4 = slice(a_prim4,CVELSTART+a_dir);
		auto bnorm4 = slice(a_prim4,CBSTART+a_dir);
		auto wnorm2 = slice(a_prim2,CVELSTART+a_dir);
		auto bnorm2 = slice(a_prim2,CBSTART+a_dir);
		// Volumetric flow rates V_u,V_b.
		auto Vu4 = Operator::_faceProduct(a_DrDetA4,wnorm4,a_DrDetA2,wnorm2,a_dir);
		auto Vb4 = Operator::_faceProduct(a_DrDetA4,bnorm4,a_DrDetA2,bnorm2,a_dir);
		auto Vu2 = a_DrDetA2*wnorm2; // Placeholder for forall.
		auto Vb2 = a_DrDetA2*bnorm2; // Placeholder for forall.

		// Advective fluxes of density, energy.
		
		auto rho4 = slice(a_prim4,CRHO);
		auto rho2 = slice(a_prim2,CRHO);
		auto fluxRho4 = Operator::_faceProduct(Vu4,rho4,Vu2,rho2,a_dir);
		auto fluxRho2 = Operator::_matrixProductAB2(Vu2,rho2);
		auto p4 = slice(a_prim4,CPRES); 
		auto p2 = slice(a_prim2,CPRES);

		// Cartesian B, velocities.
		auto w4 = slice<T,CPRIM,DIM,MEM>(a_prim4,CVELSTART);
		auto b4 = slice<T,CPRIM,DIM,MEM>(a_prim4,CBSTART);  
		auto w2 = slice<T,CPRIM,DIM,MEM>(a_prim2,CVELSTART);
		auto b2 = slice<T,CPRIM,DIM,MEM>(a_prim2,CBSTART);
		
		BoxData<T,DIM,MEM> Uface4 = Operator::_faceMatrixProductAB(a_A4,w4,a_A2,w2,a_dir);
		BoxData<T,DIM,MEM> Bface4 = Operator::_faceMatrixProductAB(a_A4,b4,a_A2,b2,a_dir);
		BoxData<T,DIM,MEM> Uface2 = Operator::_matrixProductAB2(a_A2,w2);
		BoxData<T,DIM,MEM> Bface2 = Operator::_matrixProductAB2(a_A2,b2);

		// Fluxes for velocity, magnetic fields.
		// BoxData<T,DIM,MEM> fluxuu = Operator::_faceMatrixProductAB(Uface4,fluxRho4,Uface2,fluxRho2,a_dir);
		auto fluxub = Operator::_faceMatrixProductAB(Bface4,Vu4,Bface2,Vu2,a_dir);
		// auto fluxbu = Operator::_faceMatrixProductAB(Bface4,Vb4,Uface2,Vb2,a_dir); //A mistake!
		auto fluxbu = Operator::_faceMatrixProductAB(Bface4,Vb4,Bface2,Vb2,a_dir);
		auto fluxbb = Operator::_faceMatrixProductAB(Bface4,Vb4,Bface2,Vb2,a_dir);
		fluxbb *= 1.0/(4*M_PI);


		// Fluxes for velocity, magnetic fields with Talwinder's method.
		auto wdrho4 = Operator::_faceProduct(wnorm4,rho4,wnorm2,rho2,a_dir);
		auto wdrho2 = wnorm2*rho2; // Placeholder for forall.
		BoxData<T,DIM,MEM> Uface_temp4 = Operator::_faceMatrixProductAB(a_DrDetAA4,w4,a_DrDetAA2,w2,a_dir);
		BoxData<T,DIM,MEM> Uface_temp2 = Operator::_matrixProductAB2(a_DrDetAA2,w2);
		BoxData<T,DIM,MEM> fluxuu = Operator::_faceTensorProduct(Uface_temp4,wdrho4,Uface_temp2,wdrho2,a_dir);


		
		// Energy as a function of the primitive (Cartesian) variables.
		auto Beng4 = Operator::_faceMatrixProductATB(Bface4,Bface4,Bface2,Bface2,a_dir);
		auto Beng2 = Operator::_matrixProductATB2(Bface2,Bface2);  
		Beng4 *= 1.0/(8.0*M_PI);
		Beng2 *= 1.0/(8.0*M_PI);
		
		// auto Ueng4 = Operator::_faceMatrixProductATB(Uface4,Uface4,Uface2,Uface2,a_dir);
		// auto Ueng2 = Operator::_matrixProductATB2(Uface2,Uface2);
		// Ueng4 *= .5;
		// Ueng2 *= .5; 

		//Talwinder's definition of kinetic energy
		auto w_actual4 = slice<T,CPRIM,DIM,MEM>(a_prim_actual4,CVELSTART);
		auto w_actual2 = slice<T,CPRIM,DIM,MEM>(a_prim_actual2,CVELSTART);
		auto Ueng4 = Operator::_faceMatrixProductATB(w_actual4,w_actual4,w_actual2,w_actual2,a_dir);
		auto Ueng2 = Operator::_matrixProductATB2(w_actual2,w_actual2);
		Ueng4 *= .5;
		Ueng2 *= .5;
		
		// Kinetic energy advective contribution.
		auto fluxUEng = Operator::_faceProduct(fluxRho4,Ueng4,fluxRho2,Ueng2,a_dir);
		
		// pzero = sum of thermal and magnetic pressures.
		auto pzero4 = p4 + Beng4;
		auto pzero2 = p2 + Beng2;

		// Advective + p d(1/rho) work contribution to the energy flux.
		p4 *= 1./(a_gamma - 1.0);
		p2 *= 1./(a_gamma - 1.0);
		auto thermBEng4 = Beng4 + p4 + pzero4;
		auto thermBEng2 = Beng2 + p2 + pzero2;
		
		auto fluxThermBAdv = Operator::_faceProduct(Vu4,thermBEng4,Vu2,thermBEng2,a_dir);
		
		// Non-gradient magnetic field contribution to the energy flux.
		BoxData<double,1,MEM> UDotB4 = Operator::_faceMatrixProductATB(Uface4,Bface4,Uface2,Bface2,a_dir);
		BoxData<double,1,MEM> UDotB2 = Operator::_matrixProductATB2(Uface2,Bface2);
		auto fluxBEng = Operator::_faceProduct(Vb4,UDotB4,Vb2,UDotB2,a_dir);
		fluxBEng *= (-1.0/(4.0*M_PI));

		// Pressure forces on the fluid.
		auto pForce = Operator::_faceMatrixProductAB(a_DrAdjA4,pzero4,a_DrAdjA2,pzero2,a_dir);
		
		// Assemble into flux vector.
		auto retval = forall<T,CFLUX,MEM,1>
			([ ] PROTO_LAMBDA(
							Var<T,CFLUX,MEM,1>& a_retval,
							Var<T,1,MEM>& a_fluxRho,
							Var<T,DIM,MEM>& a_fluxuu,
							Var<T,DIM,MEM>& a_fluxub,
							Var<T,DIM,MEM>& a_fluxbu,
							Var<T,DIM,MEM>& a_fluxbb,
							Var<T,1,MEM>& a_fluxThermBAdv,
							Var<T,1,MEM>& a_fluxUEng,
							Var<T,1,MEM>& a_fluxBEng,
							Var<T,DIM,MEM>& a_pforce)
			{
			a_retval(0) = a_fluxRho(0);
			for (int dir = 0; dir < DIM; dir++)
				{
				a_retval(CVELSTART+dir) = a_fluxuu(dir) + a_pforce(dir) + a_fluxbb(dir);
				a_retval(CBSTART+dir) = a_fluxub(dir) - a_fluxbu(dir);
				}
			a_retval(CENG) = a_fluxThermBAdv(0) + a_fluxUEng(0) + a_fluxBEng(0);
			},
			fluxRho4,fluxuu,fluxub,fluxbu,fluxbb,fluxThermBAdv,fluxUEng,fluxBEng,pForce);
		
		return retval;
	}



	// With Phil's operators
	/**
	 * @brief Function to calculate the right hand side of the finite volume method, Magnetic field divergence, and min dt according to CFL condition. 
	 * This function is tailored for the special mapping of the spherical grid that preserves radial flow.
	 * @param a_Rhs the output LevelBoxData.
	 * @param a_JU_ave the input LevelBoxData containing 4th order averaged product of Jacobian and conserved variables.
	 * @note If no mapping is used, J = 1.
	 */ 
	void step_spherical(LevelBoxData<double,NUMCOMPS>& a_Rhs,
			  LevelBoxData<double,NUMCOMPS>& a_JU_ave,
			  MHDLevelDataState& a_State,
			  double& a_min_dt)
	{	
		static Stencil<double> m_laplacian;
		static Stencil<double> m_deconvolve;
		static Stencil<double> m_copy;
		static Stencil<double> m_laplacian_f[DIM];
		static Stencil<double> m_deconvolve_f[DIM];
		static Stencil<double> m_convolve_f[DIM];
		static Stencil<double> m_interp_H[DIM];
		static Stencil<double> m_interp_L[DIM];
		static Stencil<double> m_interp_edge[DIM];
		static Stencil<double> m_interp_f_2nd[DIM];
		static Stencil<double> m_divergence[DIM];
		static bool initialized = false;
		if(!initialized)
		{
			m_laplacian = Stencil<double>::Laplacian();
			m_deconvolve = (-1.0/24.0)*m_laplacian + (1.0)*Shift(Point::Zeros());
			m_copy = 1.0*Shift(Point::Zeros());
			for (int dir = 0; dir < DIM; dir++)
			{
				m_laplacian_f[dir] = Stencil<double>::LaplacianFace(dir);
				m_deconvolve_f[dir] = (-1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
				m_convolve_f[dir] = (1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
				m_interp_f_2nd[dir] = 0.5*Shift(Point::Zeros()) + 0.5*Shift(-Point::Basis(dir)); 
				m_interp_H[dir] = Stencil<double>::CellToFaceH(dir);
				m_interp_L[dir] = Stencil<double>::CellToFaceL(dir);
				m_interp_edge[dir] = Stencil<double>::CellToFace(dir);
				m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
			}
			initialized =  true;
		}

		using namespace std;
		
		double a_dx = a_State.m_dx;
		double a_dy = a_State.m_dy;
		double a_dz = a_State.m_dz;
		double gamma = a_State.m_gamma;
		double dxd[3] = {a_dx, a_dy, a_dz}; // Because now its r, theta, phi
		double dt_new;
		for (auto dit : a_State.m_U){
			Box dbx0 = a_JU_ave[dit].box();
			Box dbx1 = dbx0.grow(NGHOST-NGHOST);
			a_Rhs[dit].setVal(0.0);
			Vector a_U_Sph_ave(dbx0), a_U_Sph_actual_ave(dbx0);
			MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_ave, a_JU_ave[dit], a_State.m_detAA_inv_avg[dit], a_State.m_A_inv_avg[dit], a_State.m_r2rdot_avg[dit], a_State.m_detA_avg[dit], a_State.m_A_row_mag_avg[dit], false, 4);
			MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_actual_ave, a_JU_ave[dit], a_State.m_detAA_inv_avg[dit], a_State.m_A_inv_avg[dit], a_State.m_r2rdot_avg[dit], a_State.m_detA_avg[dit], a_State.m_A_row_mag_avg[dit], false, 4);
			MHD_Mapping::Correct_V_theta_phi_at_poles(a_U_Sph_ave, a_dx, a_dy, a_dz);
			Vector W_bar = forall<double,NUMCOMPS>(consToPrimSph, a_U_Sph_ave, a_U_Sph_actual_ave, gamma);
			Vector W_ave = m_copy(W_bar);
			HDF5Handler h5;
			if (!a_State.m_min_dt_calculated){ 
				MHD_CFL::Min_dt_calc_func(dt_new, W_ave, dbx0, a_dx, a_dy, a_dz, gamma);	
				if (dt_new < a_min_dt) a_min_dt = dt_new;
			}
			
			BoxData<double,DIM,MEM,DIM> A4(dbx0);
			BoxData<double,DIM,MEM,DIM> DrDetAA4(dbx0);

			for (int d = 0; d < DIM; d++)
			{
				Vector W_ave_low_temp(dbx0), W_ave_high_temp(dbx0);
				Vector W_ave_low(dbx0), W_ave_high(dbx0);
				Vector W_ave_low_actual(dbx0), W_ave_high_actual(dbx0);

				W_ave_low_temp = m_interp_L[d](W_ave);
				W_ave_high_temp = m_interp_H[d](W_ave);

				MHD_Limiters::MHD_Limiters_4O(W_ave_low,W_ave_high,W_ave_low_temp,W_ave_high_temp,W_ave,W_bar,d,a_dx, a_dy, a_dz);			
				MHD_Mapping::W_Sph_to_W_normalized_sph(W_ave_high_actual, W_ave_low, a_State.m_A_row_mag_1_avg[dit], a_State.m_A_row_mag_2_avg[dit], a_State.m_A_row_mag_3_avg[dit], d);
				

				// Box bxFace = dbx0.grow(NGHOST).extrude(d);
				Vector F_ave_f(dbx0);
				F_ave_f.setVal(0.0);
				double dx_d = dxd[d];
				// F_ave_f = Operator::MHDSphericalFlux<double,8,8,MEMTYPE_DEFAULT>
                //   (W_ave_low,W_ave_low,A4,A4,a_State.m_Dr_detA_avg[d][dit],a_State.m_Dr_detA_avg[d][dit],a_State.m_Dr_adjA_avg[d][dit], a_State.m_Dr_adjA_avg[d][dit],gamma,d);

				if (d==0) {
					MHD_Mapping::Nineto33(A4, a_State.m_A_1_avg[dit]);
					MHD_Mapping::Nineto33(DrDetAA4, a_State.m_r2detAA_1_avg[dit]);
					F_ave_f = MHDSphericalFlux<double,8,8,MEMTYPE_DEFAULT>
                    (W_ave_high,W_ave_high,W_ave_high_actual,W_ave_high_actual,DrDetAA4,DrDetAA4,A4,A4,a_State.m_r2detA_1_avg[dit],a_State.m_r2detA_1_avg[dit],a_State.m_r2detAn_1_avg[dit], a_State.m_r2detAn_1_avg[dit],gamma,d);
				}

				if (d==1) {
					MHD_Mapping::Nineto33(A4, a_State.m_A_2_avg[dit]);
					MHD_Mapping::Nineto33(DrDetAA4, a_State.m_rrdotdetAA_2_avg[dit]);
					F_ave_f = MHDSphericalFlux<double,8,8,MEMTYPE_DEFAULT>
                    (W_ave_high,W_ave_high,W_ave_high_actual,W_ave_high_actual,DrDetAA4,DrDetAA4,A4,A4,a_State.m_rrdotdetA_2_avg[dit],a_State.m_rrdotdetA_2_avg[dit],a_State.m_rrdotd3ncn_2_avg[dit], a_State.m_rrdotd3ncn_2_avg[dit],gamma,d);
				}

				if (d==2) {
					MHD_Mapping::Nineto33(A4, a_State.m_A_3_avg[dit]);
					MHD_Mapping::Nineto33(DrDetAA4, a_State.m_rrdotdetAA_3_avg[dit]);
					F_ave_f = MHDSphericalFlux<double,8,8,MEMTYPE_DEFAULT>
                    (W_ave_high,W_ave_high,W_ave_high_actual,W_ave_high_actual,DrDetAA4,DrDetAA4,A4,A4,a_State.m_rrdotdetA_3_avg[dit],a_State.m_rrdotdetA_3_avg[dit],a_State.m_rrdotncd2n_3_avg[dit], a_State.m_rrdotncd2n_3_avg[dit],gamma,d);
				}

				Vector Rhs_d = m_divergence[d](F_ave_f);
				// cout << Rhs_d.box();
				Rhs_d *= -1./dx_d;
				a_Rhs[dit] += Rhs_d;
			}
			// Vector a_Rhs_sph(dbx0);
			// MHD_Mapping::JU_to_W_Sph_ave_calc_func(a_Rhs_sph, a_Rhs[dit], (a_State.m_detAA_inv_avg)[ dit], (a_State.m_A_inv_avg)[ dit], (a_State.m_r2rdot_avg)[ dit], (a_State.m_detA_avg)[ dit], (a_State.m_A_row_mag_avg)[ dit], inputs.gamma, true);
			// MHD_Mapping::JU_to_W_Sph_ave_calc_func2(a_Rhs_sph, a_Rhs[dit], (a_State.m_detAA_inv_avg)[ dit], (a_State.m_A_inv_avg)[ dit], (a_State.m_r2rdot_avg)[ dit], (a_State.m_detA_avg)[ dit], (a_State.m_A_row_mag_avg)[ dit], inputs.gamma, true);
			// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, a_Rhs[dit], "a_Rhs_bad");
			// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, a_Rhs_sph, "a_Rhs_sph_bad");
		}
	}






	/**
	 * @brief Function to make pressure epsilon if it goes to negative for some reason. 
	 * @param a_U the input LevelBoxData containing conserved variables.
	 * @param a_gamma gamma.
	 */ 
	PROTO_KERNEL_START
	void Fix_negative_P_calcF(const Point& a_pt,
						State&         a_U,
	            		double a_gamma)
	{
		if (inputs.Spherical_2nd_order == 1 && inputs.grid_type_global == 2){
			double rho = a_U(0);
			double v2 = 0.0;
			double B2 = 0.0;
			double gamma = a_gamma;

			for (int i = 1; i <= DIM; i++)
			{
				double v, B;
				v = a_U(i) / rho;
				B = a_U(DIM+1+i);
				v2 += v*v;
				B2 += B*B;
			}

			double p = std::max((a_U(NUMCOMPS-1-DIM) - .5 * rho * v2  - B2/8.0/c_PI) * (gamma - 1.0),1.0e-14);
			a_U(NUMCOMPS-1-DIM) = p/(gamma-1.0) + .5 * rho * v2  + B2/8.0/c_PI;
		}
	}
	PROTO_KERNEL_END(Fix_negative_P_calcF, Fix_negative_P_calc)

	/**
	 * @brief Function to make pressure epsilon if it goes to negative for some reason. 
	 * @param a_U the input LevelBoxData containing conserved variables.
	 * @param a_gamma gamma.
	 */ 
	void Fix_negative_P(BoxData<double,NUMCOMPS>& a_U,
	                    const double gamma)
	{
		forallInPlace_p(Fix_negative_P_calc, a_U, gamma);
	}

	/**
	 * @brief Function to scale flux with face area. 
	 * @param a_F_scaled the output scaled flux.
	 * @param a_F the input flux.
	 * @param a_face_area the input face area.
	 * @param a_d direction index.
	 */ 
	PROTO_KERNEL_START
	void Scale_with_A_Ff_calcF(const Point& a_pt,
						Var<double,NUMCOMPS>& a_F_scaled,
						Var<double,NUMCOMPS>& a_F,
	                   	Var<double,DIM>& a_face_area,
	                   	int a_d)
	{
		double area = a_face_area(a_d);

		for (int i = 0; i < NUMCOMPS; i++){
			a_F_scaled(i) = -a_F(i)*area;
		}
	}
	PROTO_KERNEL_END(Scale_with_A_Ff_calcF, Scale_with_A_Ff_calc)


	/**
	 * @brief Function to scale B at face with face area. 
	 * @param a_B_scaled the output scaled B_face.
	 * @param a_W_sph1 the input primitive variables on low side.
	 * @param a_W_sph2 the input primitive variables on high side.
	 * @param a_face_area the input face area.
	 * @param a_d direction index.
	 */ 
	PROTO_KERNEL_START
	void Scale_with_A_Bf_calcF(const Point& a_pt,
						Var<double,1>& a_B_scaled,
						const State& a_W_sph1,
	               		const State& a_W_sph2,
	                   	Var<double,DIM>& a_face_area,
	                   	int a_d)
	{
		double area = a_face_area(a_d);
		double B_face = 0.5*(a_W_sph1(2+DIM+a_d)+a_W_sph2(2+DIM+a_d));
		a_B_scaled(0) = -B_face*area;
		
	}
	PROTO_KERNEL_END(Scale_with_A_Bf_calcF, Scale_with_A_Bf_calc)

	/**
	 * @brief Function to scale flux at a face with cell volume. 
	 * @param a_F_scaled the output scaled flux.
	 * @param a_F the input flux.
	 * @param a_cell_volume the input cell volume.
	 */ 
	PROTO_KERNEL_START
	void Scale_with_V_calcF(const Point& a_pt,
						Var<double,NUMCOMPS>& a_F_scaled,
						Var<double,NUMCOMPS>& a_F,
	                   	Var<double,1>& a_cell_volume)
	{
		double volume = a_cell_volume(0);

		for (int i = 0; i < NUMCOMPS; i++){
			a_F_scaled(i) = a_F(i)/volume;
		}
	}
	PROTO_KERNEL_END(Scale_with_V_calcF, Scale_with_V_calc)

	/**
	 * @brief Function to scale a scalar flux at a face with cell volume. 
	 * @param a_F_scaled the output scaled flux.
	 * @param a_F the input flux.
	 * @param a_cell_volume the input cell volume.
	 */ 
	PROTO_KERNEL_START
	void Scale_with_V2_calcF(const Point& a_pt,
						Var<double,NUMCOMPS>& a_F_scaled,
						Var<double,1>& a_F,
	                   	Var<double,1>& a_cell_volume)
	{
		double volume = a_cell_volume(0);

		for (int i = 0; i < NUMCOMPS; i++){
			a_F_scaled(i) = a_F(0)/volume;
		}
	}
	PROTO_KERNEL_END(Scale_with_V2_calcF, Scale_with_V2_calc)


	// PROTO_KERNEL_START
	// void Scale_Ff_calc2F(const Point& a_pt,
	// 					Var<double,NUMCOMPS>& a_F_scaled,
	// 					Var<double,NUMCOMPS>& a_F,
	//                     Var<double,DIM>& a_x_sph_cc,
	//                    	Var<double,DIM>& a_x_sph_fc,
	//                    	Var<double,DIM>& a_dx_sp,
	//                    	int a_d)
	// {
	// 	double r_cc = a_x_sph_cc(0);
	// 	double theta_cc = a_x_sph_cc(1);

	// 	double r_fc = a_x_sph_fc(0);
	// 	double theta_fc = a_x_sph_fc(1);

	// 	double dr = a_dx_sp(0);
	// 	double dtheta = a_dx_sp(1);
	// 	double dphi = a_dx_sp(2);

	// 	if (a_d == 0){
	// 		for (int i = 0; i < NUMCOMPS; i++){
	// 			// a_F_scaled(i) = -a_F(i)*(r_fc*r_fc)/(dr*r_cc*r_cc);
	// 			a_F_scaled(i) = (r_fc*r_fc)/(dr*r_cc*r_cc);
	// 		}
	// 	}

	// 	if (a_d == 1){
	// 		for (int i = 0; i < NUMCOMPS; i++){
	// 			// a_F_scaled(i) = -a_F(i)*sin(theta_fc)/(r_cc*sin(theta_cc)*sin(dtheta));
	// 			a_F_scaled(i) = sin(theta_fc)/(r_cc*sin(theta_cc)*sin(dtheta));
	// 		}
	// 	}

	// 	if (a_d == 2){
	// 		for (int i = 0; i < NUMCOMPS; i++){
	// 			// a_F_scaled(i) = -a_F(i)/(r_cc*sin(theta_cc)*dphi);
	// 			a_F_scaled(i) = 1.0/(r_cc*sin(theta_cc)*dphi);
	// 		}
	// 	}

	// }
	// PROTO_KERNEL_END(Scale_Ff_calc2F, Scale_Ff_calc2)

	/**
	 * @brief Function to calculate Powell term for div B calculation. 
	 * @param a_P the output powell terms.
	 * @param a_W the input primitive variables.
	 * @param a_F the input divB.
	 * @param a_cell_volume the input cell volume.
	 */
	PROTO_KERNEL_START
	void Powell_Sph_2OF(const Point& a_pt,
				 State&         a_P,
	             const State&   a_W,
				 Var<double,1>& a_F,
				 Var<double,1>& a_cell_volume)
	{
	double volume = a_cell_volume(0);
#if DIM==2
		a_P(0) = (a_F(0)/volume)*0.;
		a_P(1) = (a_F(0)/volume)*a_W(4)/4.0/c_PI;
		a_P(2) = (a_F(0)/volume)*a_W(5)/4.0/c_PI;
		a_P(3) = (a_F(0)/volume)*(a_W(1)*a_W(4)/4.0/c_PI + a_W(2)*a_W(5)/4.0/c_PI);
		a_P(4) = (a_F(0)/volume)*a_W(1);
		a_P(5) = (a_F(0)/volume)*a_W(2);
#endif

#if DIM==3
		a_P(0) = (a_F(0)/volume)*0.;
		a_P(1) = (a_F(0)/volume)*a_W(5)/4.0/c_PI;
		a_P(2) = (a_F(0)/volume)*a_W(6)/4.0/c_PI;
		a_P(3) = (a_F(0)/volume)*a_W(7)/4.0/c_PI;
		a_P(4) = (a_F(0)/volume)*(a_W(1)*a_W(5)/4.0/c_PI + a_W(2)*a_W(6)/4.0/c_PI + a_W(3)*a_W(7)/4.0/c_PI);
		a_P(5) = (a_F(0)/volume)*a_W(1);
		a_P(6) = (a_F(0)/volume)*a_W(2);
		a_P(7) = (a_F(0)/volume)*a_W(3);

#endif
	}
	PROTO_KERNEL_END(Powell_Sph_2OF, Powell_Sph_2O)



	/**
	 * @brief Function to calculate B at face. 
	 * @param a_B_face the output B a face.
	 * @param a_W_sph1 the input primitive variables on low side.
	 * @param a_W_sph2 the input primitive variables on high side.
	 * @param a_d the input direction.
	 */
	PROTO_KERNEL_START
	void B_face_calcF(const Point& a_pt,
                    Var<double,1>& a_B_face,
	               const State& a_W_sph1,
	               const State& a_W_sph2,
	               int a_d)
	{
		a_B_face(0) = 0.5*(a_W_sph1(2+DIM+a_d)+a_W_sph2(2+DIM+a_d));
	}
	PROTO_KERNEL_END(B_face_calcF, B_face_calc)


	/**
	 * @brief Function to calculate the right hand side of the finite volume method, Magnetic field divergence, and min dt according to CFL condition. 
	 * @param a_Rhs the output LevelBoxData.
	 * @param a_U the input LevelBoxData containing conserved variables.
	 */
	void step_spherical_2O(LevelBoxData<double,NUMCOMPS>& a_Rhs,
			  LevelBoxData<double,NUMCOMPS>& a_U,
			  MHDLevelDataState& a_State,
			  double& a_min_dt)
	{
		
		using namespace std;
		double a_dx = a_State.m_dx;
		double a_dy = a_State.m_dy;
		double a_dz = a_State.m_dz;
		double gamma = a_State.m_gamma;
		double dxd[3] = {a_dx, a_dy, a_dz};
		
		static Stencil<double> m_divergence[DIM];
		static bool initialized = false;
		if(!initialized)
		{
			for (int dir = 0; dir < DIM; dir++)
			{
				m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
			}
			initialized =  true;
		}
		double dt_new;
		for (auto dit : a_State.m_U){	
			Box dbx0 = a_U[dit].box();
			Box dbx1 = a_State.m_U[dit].box();
			Vector F_f_sph(dbx0), F_f(dbx0), F_f_scaled(dbx0), Rhs_d(dbx0), RhsV(dbx0);
			Scalar B_f_sph(dbx0), B_f_scaled(dbx0), Rhs_d_divB(dbx0), RhsV_divB(dbx0);
			RhsV.setVal(0.0);
			RhsV_divB.setVal(0.0);

			MHDOp::Fix_negative_P(a_U[dit],inputs.gamma);	
			Vector W_low_temp(dbx0), W_high_temp(dbx0), W_low(dbx0), W_high(dbx0);
			BoxData<double,NUMCOMPS> W_sph(dbx0);
			Vector W_cart  = forall<double,NUMCOMPS>(consToPrim, a_U[dit], gamma);
			MHD_Mapping::Cartesian_to_Spherical(W_sph, W_cart, a_State.m_x_sph_cc[dit]);
			MHD_Mapping::Correct_V_theta_phi_at_poles(W_sph, a_dx, a_dy, a_dz);	

			if (!a_State.m_min_dt_calculated){ 
				MHD_CFL::Min_dt_calc_func(dt_new, W_sph, dbx1, a_dx, a_dy, a_dz, gamma);	
				if (dt_new < a_min_dt) a_min_dt = dt_new;
			}

			for (int d = 0; d < DIM; d++)
			{
				MHD_Limiters::MHD_Limiters_minmod(W_low,W_high,W_sph,a_State.m_x_sph_cc[dit],a_State.m_dx_sph[dit],d);
				MHD_Riemann_Solvers::Roe8Wave_Solver(F_f_sph,W_low,W_high,d,gamma);
				if (d==0) MHD_Mapping::Spherical_to_Cartesian(F_f, F_f_sph, a_State.m_x_sph_fc_1[dit]);
				if (d==1) MHD_Mapping::Spherical_to_Cartesian(F_f, F_f_sph, a_State.m_x_sph_fc_2[dit]);
				if (d==2) MHD_Mapping::Spherical_to_Cartesian(F_f, F_f_sph, a_State.m_x_sph_fc_3[dit]);
				forallInPlace_p(Scale_with_A_Ff_calc, F_f_scaled, F_f, a_State.m_face_area[dit], d);
				Rhs_d = m_divergence[d](F_f_scaled);
				RhsV += Rhs_d;

				if (!a_State.m_divB_calculated){
					forallInPlace_p(Scale_with_A_Bf_calc, B_f_scaled, W_low,W_high, a_State.m_face_area[dit], d);
					Rhs_d_divB = m_divergence[d](B_f_scaled);
					RhsV_divB += Rhs_d_divB;
				}
			}
			forallInPlace_p(Scale_with_V_calc, a_Rhs[dit], RhsV, a_State.m_cell_volume[dit]);
			if (!a_State.m_divB_calculated) forallInPlace_p(Powell_Sph_2O,a_State.m_divB[dit],W_cart,RhsV_divB,a_State.m_cell_volume[dit]);		
		}
	}
}
