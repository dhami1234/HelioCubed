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
#include "MHD_Operator.H"

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
			if (!a_State.m_divB_calculated && inputs.takedivBstep == 1) a_State.m_divB[dit].setVal(0.0);

			HDF5Handler h5;

			auto a_U_ave = Operator::_cellTensorQuotient(a_JU_ave[dit],a_State.m_J[dit],a_JU_ave[dit],a_State.m_J[dit]);

			Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U_ave, gamma);
			Vector U = m_deconvolve(a_U_ave);
			Vector W  = forall<double,NUMCOMPS>(consToPrim,U, gamma);
			Vector W_ave = m_laplacian(W_bar,1.0/24.0);
			W_ave += W;
			if (!a_State.m_min_dt_calculated){ 
				MHD_CFL::Min_dt_calc_func(dt_new, W_ave, dbx2, a_dx, a_dy, a_dz, gamma);	
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
				if (!a_State.m_divB_calculated && inputs.takedivBstep == 1) {BF_f_mapped.setVal(0.0);}
				double dx_d = dxd[d];
				for (int s = 0; s < DIM; s++) {
					if (inputs.Riemann_solver_type == 1) {
						MHD_Riemann_Solvers::Rusanov_Solver(F_f_temp,W_low,W_high,s,gamma);
					}
					if (inputs.Riemann_solver_type == 2) {
						MHD_Riemann_Solvers::Roe8Wave_Solver(F_f_temp,W_low,W_high,s,gamma);
					}
					forallInPlace_p(Fill_flux_calc, F_f, F_f_temp, s);
					if (!a_State.m_divB_calculated && inputs.takedivBstep == 1){
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

				if (!a_State.m_divB_calculated && inputs.takedivBstep == 1){
					BoxData<double,1,MEM,NUMCOMPS> Bfluxdir = Operator::_faceMatrixProductATB(a_State.m_NT[d][dit],BF_ave_f,a_State.m_NT[d][dit],BF_ave_f,d);
					forallInPlace_p(Transpose_calc, BF_f_mapped, Bfluxdir);
					Vector B_Rhs_d = m_divergence[d](BF_f_mapped);
					// B_Rhs_d *= -1./pow(dx_d,DIM);
					B_Rhs_d *= -1./volume;
					a_State.m_divB[dit] += B_Rhs_d;
				}
			}
			if (!a_State.m_divB_calculated && inputs.takedivBstep == 1){
				Vector Powell_term = forall<double,NUMCOMPS>(Powell,W_ave);
				a_State.m_divB[dit] *= Powell_term;
			}
		}
	}

	PROTO_KERNEL_START
	void QuotientF(State&         a_RHSbyJ,
	             const State&   a_RHS,
				 const Var<double,1>& a_J)
	{
		for (int i=0;i<NUMCOMPS;i++){
			a_RHSbyJ(i) = a_RHS(i)/a_J(0);
		}
	}
	PROTO_KERNEL_END(QuotientF, Quotient)

	PROTO_KERNEL_START
	void ProductF(State&         a_a,
	             const State&   a_b,
				 const Var<double,1>& a_c)
	{
		for (int i=0;i<NUMCOMPS;i++){
			a_a(i) = a_b(i)*a_c(0);
		}
	}
	PROTO_KERNEL_END(ProductF, Product)

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
			HDF5Handler h5;
			Box dbx0 = a_JU_ave[dit].box();
			Box dbx1 = dbx0.grow(0-NGHOST);
			a_Rhs[dit].setVal(0.0);
			Vector a_U_Sph_ave(dbx0), a_U_Sph_actual_ave(dbx0), a_U_cart_ave(dbx0), a_U_cart_ave_temp(dbx0), Powell_term(dbx0), F_ave_f(dbx0);
			Scalar RhsV_divB(dbx0), NTB(dbx0), NTV(dbx0);
			if (!a_State.m_divB_calculated && inputs.takedivBstep == 1) RhsV_divB.setVal(0.0);
			if (!a_State.m_divV_calculated && inputs.non_linear_visc_apply == 1) a_State.m_divV[dit].setVal(0.0);
			MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_ave, a_JU_ave[dit], a_State.m_detAA_inv_avg[dit], a_State.m_r2rdot_avg[dit], a_State.m_detA_avg[dit], a_State.m_A_row_mag_avg[dit], false, 2);
			MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_actual_ave, a_JU_ave[dit], a_State.m_detAA_inv_avg[dit], a_State.m_r2rdot_avg[dit], a_State.m_detA_avg[dit], a_State.m_A_row_mag_avg[dit], true, 2);
			MHD_Mapping::JU_to_U_ave_calc_func(a_U_cart_ave, a_JU_ave[dit], a_State.m_r2rdot_avg[dit], a_State.m_detA_avg[dit]);
			MHD_Mapping::Correct_V_theta_phi_at_poles(a_U_Sph_ave, a_dx, a_dy, a_dz);
			MHD_Mapping::Correct_V_theta_phi_at_poles(a_U_Sph_actual_ave, a_dx, a_dy, a_dz);
			Vector W_bar = forall<double,NUMCOMPS>(consToPrimSph, a_U_Sph_ave, a_U_Sph_actual_ave, gamma);
			Vector W_bar_actual = forall<double,NUMCOMPS>(consToPrimSph, a_U_Sph_actual_ave, a_U_Sph_actual_ave, gamma);
			
			Vector U = m_deconvolve(a_U_Sph_ave);
			Vector U_actual = m_deconvolve(a_U_Sph_actual_ave);
			Vector W = forall<double,NUMCOMPS>(consToPrimSph, U, U_actual, gamma);
			Vector W_ave = m_laplacian(W_bar,1.0/24.0);
			W_ave += W;
			Vector W_actual = forall<double,NUMCOMPS>(consToPrimSph, a_U_Sph_actual_ave, a_U_Sph_actual_ave, gamma);

			
			if (!a_State.m_min_dt_calculated){ 
				// MHD_CFL::Min_dt_calc_func(dt_new, W_actual, dbx1, a_dx, a_dy, a_dz, gamma);	
				MHD_CFL::Min_dt_calc_func(dt_new, W_bar_actual, dbx1, a_dx, a_dy, a_dz, gamma);	
				if (dt_new < a_min_dt) a_min_dt = dt_new;
			}

			for (int d = 0; d < DIM; d++)
			{
				Vector W_ave_low_temp(dbx0), W_ave_high_temp(dbx0);
				Vector W_ave_low(dbx0), W_ave_high(dbx0);
				Vector W_ave_low_actual(dbx0), W_ave_high_actual(dbx0);

				// W_ave_low_temp = m_interp_L[d](W_ave);
				// W_ave_high_temp = m_interp_H[d](W_ave);
				// MHD_Limiters::MHD_Limiters_4O(W_ave_low,W_ave_high,W_ave_low_temp,W_ave_high_temp,W_ave,W_bar,d,a_dx, a_dy, a_dz);		
				
				MHD_Limiters::MHD_Limiters_minmod(W_ave_low,W_ave_high,W_bar,a_State.m_x_sph_cc[dit],a_State.m_dx_sph[dit],d);
				
				MHD_Mapping::W_Sph_to_W_normalized_sph(W_ave_low_actual, W_ave_low, a_State.m_A_row_mag_face_avg[d][dit], d);
				MHD_Mapping::W_Sph_to_W_normalized_sph(W_ave_high_actual, W_ave_high, a_State.m_A_row_mag_face_avg[d][dit], d);
				// Vector W_low_boxed(dbx1);
				// W_ave_low_actual.copyTo(W_low_boxed);
				// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, W_low_boxed, "W_low_4O_"+to_string(d));
				// F_ave_f.setVal(0.0);
				auto bnorm4 = slice(W_ave_low,CBSTART+d)+slice(W_ave_high,CBSTART+d);
				bnorm4 *= 0.5;
				auto vnorm4 = slice(W_ave_low,CVELSTART+d)+slice(W_ave_high,CVELSTART+d);
				vnorm4 *= 0.5;
				double dx_d = dxd[d];	
				MHD_Riemann_Solvers::Spherical_Riemann_Solver(F_ave_f, W_ave_low, W_ave_high, W_ave_low_actual, W_ave_high_actual, a_State.m_Dr_detA_avg[d][dit], a_State.m_Dr_detA_A_avg[d][dit], a_State.m_Dr_AdjA_avg[d][dit], a_State.m_A_row_mag_face_avg[d][dit], d, gamma, a_dx, a_dy, a_dz, 4);	
				// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, F_ave_f, "F_ave_f_"+to_string(d));
				Vector Rhs_d = m_divergence[d](F_ave_f);
				Rhs_d *= -1./dx_d;
				a_Rhs[dit] += Rhs_d;
				if (!a_State.m_divB_calculated && inputs.takedivBstep == 1){
					// NTB = Operator::_faceProduct(a_State.m_Dr_detA_avg[d][dit],bnorm4,a_State.m_Dr_detA_avg[d][dit],bnorm4,d);
					NTB = a_State.m_Dr_detA_avg[d][dit]*bnorm4;
					Scalar B_Rhs_d = m_divergence[d](NTB);
					B_Rhs_d *= -1./dx_d;
					RhsV_divB += B_Rhs_d;
				}
				if (!a_State.m_divV_calculated && inputs.non_linear_visc_apply == 1){
					// NTV = Operator::_faceProduct(a_State.m_Dr_detA_avg[d][dit],vnorm4,a_State.m_Dr_detA_avg[d][dit],vnorm4,d);
					NTV = a_State.m_Dr_detA_avg[d][dit]*vnorm4;
					Scalar V_Rhs_d = m_divergence[d](NTV);
					V_Rhs_d *= -1./dx_d;
					a_State.m_divV[dit] += V_Rhs_d;
				}
			}

			if (!a_State.m_divB_calculated && inputs.takedivBstep == 1){
				
				Vector W_cart = forall<double,NUMCOMPS>(consToPrim,a_U_cart_ave, gamma);
				Powell_term = forall<double,NUMCOMPS>(Powell,W_cart);
				a_State.m_divB[dit] = forall<double,NUMCOMPS>(Product,Powell_term,RhsV_divB);
				// a_State.m_divB[dit] = MHD_Operator::_cellTensorProduct(Powell_term,RhsV_divB,Powell_term,RhsV_divB);
			}
			// Vector RHSbyJ = forall<double,NUMCOMPS>(Quotient,a_State.m_divB[dit],a_State.m_Jacobian_ave[dit]);
			// Vector RHSbyJ = forall<double,NUMCOMPS>(Quotient,a_Rhs[dit],a_State.m_Jacobian_ave[dit]);
			// Vector RHS(dbx1);
			// RHSbyJ.copyTo(RHS);
			// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, RHS, "RHS_4");
		}
	}

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
			// Box dbx1 = a_State.m_U[dit].box();
			Box dbx1 = dbx0.grow(0-NGHOST);
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
			HDF5Handler h5;

			// Vector W_boxed(dbx1);
			// W_sph.copyTo(W_boxed);
			// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, W_boxed, "W_2");

			if (!a_State.m_min_dt_calculated){ 
				MHD_CFL::Min_dt_calc_func(dt_new, W_sph, dbx1, a_dx, a_dy, a_dz, gamma);	
				if (dt_new < a_min_dt) a_min_dt = dt_new;
			}
			// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, W_sph, "W_ave2O");
			for (int d = 0; d < DIM; d++)
			{
				MHD_Limiters::MHD_Limiters_minmod(W_low,W_high,W_sph,a_State.m_x_sph_cc[dit],a_State.m_dx_sph[dit],d);
				// Vector W_low_boxed(dbx1);
				// W_low.copyTo(W_low_boxed);
				// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, W_low_boxed, "W_low_2O_"+to_string(d));
				if (inputs.Riemann_solver_type == 1) MHD_Riemann_Solvers::Rusanov_Solver(F_f_sph,W_low,W_high,d,gamma);
				if (inputs.Riemann_solver_type == 2) MHD_Riemann_Solvers::Roe8Wave_Solver(F_f_sph,W_low,W_high,d,gamma);
				if (d==0) MHD_Mapping::Spherical_to_Cartesian(F_f, F_f_sph, a_State.m_x_sph_fc_1[dit]);
				if (d==1) MHD_Mapping::Spherical_to_Cartesian(F_f, F_f_sph, a_State.m_x_sph_fc_2[dit]);
				if (d==2) MHD_Mapping::Spherical_to_Cartesian(F_f, F_f_sph, a_State.m_x_sph_fc_3[dit]);
				// Vector F_f_boxed(dbx1);
				// F_f.copyTo(F_f_boxed);
				// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, F_f_boxed, "F_f_2O_"+to_string(d));
				forallInPlace_p(Scale_with_A_Ff_calc, F_f_scaled, F_f, a_State.m_face_area[dit], d);
				Rhs_d = m_divergence[d](F_f_scaled);
				RhsV += Rhs_d;

				if (!a_State.m_divB_calculated && inputs.takedivBstep == 1){
					forallInPlace_p(Scale_with_A_Bf_calc, B_f_scaled, W_low,W_high, a_State.m_face_area[dit], d);
					Rhs_d_divB = m_divergence[d](B_f_scaled);
					RhsV_divB += Rhs_d_divB;
				}
			}
			forallInPlace_p(Scale_with_V_calc, a_Rhs[dit], RhsV, a_State.m_cell_volume[dit]);
			if (!a_State.m_divB_calculated && inputs.takedivBstep == 1) forallInPlace_p(Powell_Sph_2O,a_State.m_divB[dit],W_cart,RhsV_divB,a_State.m_cell_volume[dit]);		
			// Vector RHS(dbx1);
			// a_State.m_divB[dit].copyTo(RHS);
			// a_Rhs[dit].copyTo(RHS);
			// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, RHS, "RHS_2");
		}
	}
}
