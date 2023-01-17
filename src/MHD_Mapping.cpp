#include "Proto.H"
#include "MHD_Mapping.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
// #include "Proto_WriteBoxData.H"
#include "MHDOp.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Constants.H"
//#include "MHDLevelDataRK4.H"
extern Parsefrominputs inputs;

// constexpr MemType MEM = MEMTYPE_DEFAULT;

typedef BoxData<double,1,HOST> Scalar;
typedef BoxData<double,NUMCOMPS,HOST> Vector;
/// @brief MHD_Mapping namespace
namespace MHD_Mapping {

	PROTO_KERNEL_START
	void iotaFuncF(Point           & a_p,
	               V               & a_X,
	               const double a_dx,
	               const double a_dy,
	               const double a_dz)
	{
		for (int ii = 0; ii < DIM; ii++)
		{
			double dxd;
			if (ii == 0) dxd = a_dx;
			if (ii == 1) dxd = a_dy;
			if (ii == 2) dxd = a_dz;
			a_X(ii) = a_p[ii]*dxd + 0.5*dxd;
		}
	}
	PROTO_KERNEL_END(iotaFuncF,iotaFunc)

	PROTO_KERNEL_START
	void iotaFuncCornerF(Point           & a_p,
	               V               & a_X,
	               const double a_dx,
	               const double a_dy,
	               const double a_dz)
	{
		for (int ii = 0; ii < DIM; ii++)
		{
			double dxd;
			if (ii == 0) dxd = a_dx;
			if (ii == 1) dxd = a_dy;
			if (ii == 2) dxd = a_dz;
			a_X(ii) = a_p[ii]*dxd;
		}
	}
	PROTO_KERNEL_END(iotaFuncCornerF,iotaFuncCorner)

	PROTO_KERNEL_START
	void iotaFuncFaceF(Point           & a_p,
	                   V               & a_X,
	                   const double a_dx,
	                   const double a_dy,
	                   const double a_dz,
	                   int a_d)
	{
		for (int ii = 0; ii < DIM; ii++)
		{
			double dxd;
			if (ii == 0) dxd = a_dx;
			if (ii == 1) dxd = a_dy;
			if (ii == 2) dxd = a_dz;
			if (ii == a_d) {
				a_X(ii) = a_p[ii]*dxd;
			} else {
				a_X(ii) = a_p[ii]*dxd + 0.5*dxd;
			}
		}
	}
	PROTO_KERNEL_END(iotaFuncFaceF,iotaFuncFace)

	void eta_calc(BoxData<double,DIM>& a_eta,
	              const Box& a_bx,
	              const double a_dx,
	              const double a_dy,
	              const double a_dz)
	{
		forallInPlace_p(iotaFunc, a_bx, a_eta, a_dx, a_dy, a_dz);
	}

	void etaCorner_calc(BoxData<double,DIM>& a_eta,
	              const Box& a_bx,
	              const double a_dx,
	              const double a_dy,
	              const double a_dz)
	{
		forallInPlace_p(iotaFuncCorner, a_bx, a_eta, a_dx, a_dy, a_dz);
	}

	void etaFace_calc(BoxData<double,DIM>& a_eta,
	                  const Box& a_bx,
	                  const double a_dx,
	                  const double a_dy,
	                  const double a_dz,
	                  int a_d)
	{
		forallInPlace_p(iotaFuncFace, a_bx, a_eta, a_dx, a_dy, a_dz, a_d);
	}

	PROTO_KERNEL_START
	void eta_to_xF(const Point& a_pt,
					Var<double,DIM>& a_x,
	               const Var<double,DIM>& a_eta)
	{

#if DIM == 2
		if (inputs.grid_type_global == 0) {
			a_x(0) = a_eta(0);
			a_x(1) = a_eta(1);
		}
		if (inputs.grid_type_global == 1) {
			double C1 = inputs.C1_fix;
			double C2 = C1;
			a_x(0) = a_eta(0) + C1*sin(2.0*c_PI*a_eta(0))*sin(2.0*c_PI*a_eta(1));
			a_x(1) = a_eta(1) + C2*sin(2.0*c_PI*a_eta(0))*sin(2.0*c_PI*a_eta(1));
		}
		if (inputs.grid_type_global == 2) {
			a_x(0) = (inputs.r_in*c_AU + a_eta(0)*(inputs.r_out*c_AU-inputs.r_in*c_AU))*cos(2.0*c_PI*a_eta(1));
			a_x(1) = (inputs.r_in*c_AU + a_eta(0)*(inputs.r_out*c_AU-inputs.r_in*c_AU))*sin(2.0*c_PI*a_eta(1));
		}
#endif
#if DIM == 3
		if (inputs.grid_type_global == 0) {
			a_x(0) = a_eta(0);
			a_x(1) = a_eta(1);
			a_x(2) = a_eta(2);
		}
		if (inputs.grid_type_global == 2) {
			if (inputs.Spherical_2nd_order == 0){
				double R_t = (inputs.r_out*c_AU - inputs.r_in*c_AU)/(exp(inputs.C_rad) - 1.0);
				double r = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*a_eta(0)) - 1.0);
				a_x(0) = r*sin(c_PI*a_eta(1))*cos(2.0*c_PI*a_eta(2));
				a_x(1) = r*sin(c_PI*a_eta(1))*sin(2.0*c_PI*a_eta(2));
				a_x(2) = r*cos(c_PI*a_eta(1));
			}
			if (inputs.Spherical_2nd_order == 1){
				double R_t = (inputs.r_out*c_AU - inputs.r_in*c_AU)/(exp(inputs.C_rad) - 1.0);
				double r = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*a_eta(0)) - 1.0);
				double theta = c_PI*a_eta(1);
				double phi = 2.0*c_PI*a_eta(2);

				if (theta < 0){
					theta = -theta;
					if (phi < c_PI){
						phi = phi + c_PI;
					} else {
						phi = phi - c_PI;
					}
				}
				if (theta > c_PI){
					theta = c_PI - (theta - c_PI);
					if (phi < c_PI){
						phi = phi + c_PI;
					} else {
						phi = phi - c_PI;
					}
				}

				a_x(0) = r*sin(theta)*cos(phi);
				a_x(1) = r*sin(theta)*sin(phi);
				a_x(2) = r*cos(theta);
			}

		}

#endif
	}
	PROTO_KERNEL_END(eta_to_xF, eta_to_x)

	
	void eta_to_x_calc(BoxData<double,DIM>& a_x,
	                   const BoxData<double,DIM>& a_eta,
					   const Box& dbx0)
	{
		forallInPlace_p(eta_to_x,dbx0,a_x,a_eta);
	}

	void phys_coords_calc(BoxData<double,DIM>& a_x,
	                      const Box& dbx1,
	                      const double a_dx,
	                      const double a_dy,
	                      const double a_dz)
	{
		Box dbx0=a_x.box();
		BoxData<double,DIM> eta(dbx0);
		MHD_Mapping::eta_calc(eta,dbx0,a_dx, a_dy, a_dz);
		MHD_Mapping::eta_to_x_calc(a_x,eta, dbx0);		
	}

	void phys_coords_face_calc(BoxData<double,DIM>& a_x,
	                      const Box& dbx1,
	                      const double a_dx,
	                      const double a_dy,
	                      const double a_dz,
						  const int a_d)
	{
		Box dbx0=a_x.box();
		BoxData<double,DIM> eta(dbx0);
		MHD_Mapping::etaFace_calc(eta,dbx0,a_dx, a_dy, a_dz, a_d);
		MHD_Mapping::eta_to_x_calc(a_x,eta, dbx0);		
	}
	

	PROTO_KERNEL_START
	void a_U_demapped_2nd_order_calcF(  State& a_U_demapped,
	                                    const Var<double,1>& a_Jacobian_ave,
	                                    const State& a_a_U)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			a_U_demapped(i) = a_a_U(i)/a_Jacobian_ave(0);
		}

	}
	PROTO_KERNEL_END(a_U_demapped_2nd_order_calcF, a_U_demapped_2nd_order_calc)


	void JU_to_U_2ndOrdercalc(BoxData<double,NUMCOMPS>& a_U_demapped,
	                          const BoxData<double,NUMCOMPS>& a_a_U,
	                          const BoxData<double,1>& a_Jacobian_ave,
	                          const Box& a_dbx0)
	{
		a_U_demapped = forall<double,NUMCOMPS>(a_U_demapped_2nd_order_calc,a_Jacobian_ave,a_a_U);
	}


	void JU_to_W_bar_calc(BoxData<double,NUMCOMPS>& a_W_bar,
	                      const BoxData<double,NUMCOMPS>& a_JU,
	                      BoxData<double,1>& a_J,
	                      const double a_gamma)
	{
		Box dbx1 = a_JU.box();
		Vector a_U(dbx1);
		MHD_Mapping::JU_to_U_2ndOrdercalc(a_U, a_JU, a_J, dbx1);
		MHDOp::consToPrimcalc(a_W_bar,a_U,a_gamma);
	}



	PROTO_KERNEL_START
	void out_data_joinF(Var<double,NUMCOMPS+DIM>& a_out_data,
	                    const Var<double,DIM>& a_phys_coords,
	                    const Var<double,NUMCOMPS>& a_W)
	{
		for (int i=0; i< DIM; i++) {
			if (inputs.grid_type_global == 2){
				a_out_data(i) = a_phys_coords(i)/c_AU; // AU
			} else {
				a_out_data(i) = a_phys_coords(i); // cm
			}
		}
		
		if (inputs.grid_type_global == 2){
			a_out_data(DIM+0) = a_W(0)/c_MP; // /cm^3
			a_out_data(DIM+1) = a_W(1)/1e5; // km/s
			a_out_data(DIM+2) = a_W(2)/1e5; // km/s
			a_out_data(DIM+3) = a_W(3)/1e5; // km/s
			a_out_data(DIM+4) = a_W(4); // dyne/cm*2
			a_out_data(DIM+5) = a_W(5); // Gauss
			a_out_data(DIM+6) = a_W(6); // Gauss
			a_out_data(DIM+7) = a_W(7); // Gauss
		} else {
			a_out_data(DIM+0) = a_W(0); // g/cm^3
			a_out_data(DIM+1) = a_W(1); // cm/s
			a_out_data(DIM+2) = a_W(2); // cm/s
			a_out_data(DIM+3) = a_W(3); // cm/s
			a_out_data(DIM+4) = a_W(4); // dyne/cm*2
			a_out_data(DIM+5) = a_W(5); // Gauss
			a_out_data(DIM+6) = a_W(6); // Gauss
			a_out_data(DIM+7) = a_W(7); // Gauss
		}

	}
	PROTO_KERNEL_END(out_data_joinF, out_data_join)


	void out_data_calc(BoxData<double,NUMCOMPS+DIM>& a_out_data,
	                   const BoxData<double,DIM>& a_phys_coords,
	                   const BoxData<double,NUMCOMPS>& a_W)
	{
		a_out_data = forall<double,NUMCOMPS+DIM>(out_data_join, a_phys_coords, a_W);
	}


	PROTO_KERNEL_START
	void Spherical_map_calcF( Point& a_pt,
	                          Var<double,1>& a_Jacobian_ave,
							  Var<double,DIM*DIM>& a_A_avg,
							  Var<double,DIM*DIM>& a_A_inv_avg,
	                          Var<double,DIM*DIM>& a_A_1_avg,
	                          Var<double,DIM*DIM>& a_A_2_avg,
	                          Var<double,DIM*DIM>& a_A_3_avg,
							  Var<double,DIM*DIM>& a_A_inv_1_avg,
	                          Var<double,DIM*DIM>& a_A_inv_2_avg,
	                          Var<double,DIM*DIM>& a_A_inv_3_avg,
	                          Var<double,DIM*DIM>& a_detAA_avg,
	                          Var<double,DIM*DIM>& a_detAA_inv_avg,
	                          Var<double,1>& a_r2rdot_avg,
	                          Var<double,1>& a_detA_avg,
							  Var<double,DIM>& a_A_row_mag_avg,
	                          Var<double,1>& a_r2detA_1_avg,
	                          Var<double,DIM*DIM>& a_r2detAA_1_avg,
	                          Var<double,DIM>& a_r2detAn_1_avg,
	                          Var<double,DIM>& a_n_1_avg,
							  Var<double,DIM>& a_A_row_mag_1_avg,
	                          Var<double,1>& a_rrdotdetA_2_avg,
	                          Var<double,DIM*DIM>& a_rrdotdetAA_2_avg,
	                          Var<double,DIM>& a_rrdotd3ncn_2_avg,
							  Var<double,DIM>& a_A_row_mag_2_avg,
	                          Var<double,1>& a_rrdotdetA_3_avg,
	                          Var<double,DIM*DIM>& a_rrdotdetAA_3_avg,
	                          Var<double,DIM>& a_rrdotncd2n_3_avg,
							  Var<double,DIM>& a_A_row_mag_3_avg,
	                          const double a_dx,
	                          const double a_dy,
	                          const double a_dz,
							  bool a_exchanged_yet,
							  bool a_r_dir_turn)
	{

		double R0 = inputs.r_in*c_AU;
		double R1 = inputs.r_out*c_AU;
		double c = inputs.C_rad;
		double Rt = (R1 - R0)/(exp(c) - 1.0);

		double dE1 = a_dx;
		double dE2 = a_dy;
		double dE3 = a_dz;

		double E1, E2, E3;
/////  CELL AVERAGED VALUES

		//For the special mapping that preserves radial flow, Phil has suggested to use E2 as theta and E3 as Phi.
		
		E1 = (a_pt[0] + 0.5)*dE1;
		E2 = (a_pt[1] + 0.5)*dE2;
		E3 = (a_pt[2] + 0.5)*dE3;
		if (a_exchanged_yet){
			if (E2 > 0.0 && E2 < 1.0) return;
		}
		
		if (E2 < 0.0){
			E2 = -E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}

		if (E2 > 1.0){
			E2 = 2.0 - E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}
		
		if (!a_r_dir_turn){
			a_detAA_avg(0)	=	-0.25*((-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI))*(sin((dE3 - 2*E3)*c_PI) + sin((dE3 + 2*E3)*c_PI)))/(dE2*dE3);
			a_detAA_avg(1)	=	(c_PI*sin(dE2*c_PI)*sin(2*E2*c_PI)*(sin((dE3 - 2*E3)*c_PI) + sin((dE3 + 2*E3)*c_PI)))/(2.*dE2*dE3);
			a_detAA_avg(2)	=	(c_PI*sin(dE3*c_PI)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI))*sin(2*E3*c_PI))/(dE2*dE3);
			a_detAA_avg(3)	=	-0.5*(sin(dE3*c_PI)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI))*sin(2*E3*c_PI))/(dE2*dE3);
			a_detAA_avg(4)	=	(c_PI*sin(dE2*c_PI)*sin(dE3*c_PI)*sin(2*E2*c_PI)*sin(2*E3*c_PI))/(dE2*dE3);
			a_detAA_avg(5)	=	(c_PI*(4*dE2*c_PI*cos(2*E3*c_PI)*sin(dE3*c_PI) - (sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI))*(sin((dE3 - 2*E3)*c_PI) + sin((dE3 + 2*E3)*c_PI))))/(2.*dE2*dE3);
			a_detAA_avg(6)	=	(c_PI*sin(dE2*c_PI)*sin(2*E2*c_PI))/dE2;
			a_detAA_avg(7)	=	(pow(c_PI,2)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI)))/(2.*dE2);
			a_detAA_avg(8)	=	0;

			double det_detAA_avg = a_detAA_avg(0)*(a_detAA_avg(4)*a_detAA_avg(8) - a_detAA_avg(7)*a_detAA_avg(5)) - a_detAA_avg(1)*(a_detAA_avg(3)*a_detAA_avg(8) - a_detAA_avg(5)*a_detAA_avg(6)) + a_detAA_avg(2)*(a_detAA_avg(3)*a_detAA_avg(7) - a_detAA_avg(4)*a_detAA_avg(6));

			a_detAA_inv_avg(0) = (a_detAA_avg(4)*a_detAA_avg(8) - a_detAA_avg(5)*a_detAA_avg(7))/det_detAA_avg;
			a_detAA_inv_avg(1) = (a_detAA_avg(2)*a_detAA_avg(7) - a_detAA_avg(1)*a_detAA_avg(8))/det_detAA_avg;
			a_detAA_inv_avg(2) = (a_detAA_avg(1)*a_detAA_avg(5) - a_detAA_avg(2)*a_detAA_avg(4))/det_detAA_avg;
			a_detAA_inv_avg(3) = (a_detAA_avg(5)*a_detAA_avg(6) - a_detAA_avg(3)*a_detAA_avg(8))/det_detAA_avg;
			a_detAA_inv_avg(4) = (a_detAA_avg(0)*a_detAA_avg(8) - a_detAA_avg(2)*a_detAA_avg(6))/det_detAA_avg;
			a_detAA_inv_avg(5) = (a_detAA_avg(2)*a_detAA_avg(3) - a_detAA_avg(0)*a_detAA_avg(5))/det_detAA_avg;
			a_detAA_inv_avg(6) = (a_detAA_avg(3)*a_detAA_avg(7) - a_detAA_avg(4)*a_detAA_avg(6))/det_detAA_avg;
			a_detAA_inv_avg(7) = (a_detAA_avg(1)*a_detAA_avg(6) - a_detAA_avg(0)*a_detAA_avg(7))/det_detAA_avg;
			a_detAA_inv_avg(8) = (a_detAA_avg(0)*a_detAA_avg(4) - a_detAA_avg(1)*a_detAA_avg(3))/det_detAA_avg;


			a_A_avg(0) = (2*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI))/(dE2*dE3*pow(c_PI,2));
			a_A_avg(1) = (2*cos(E2*c_PI)*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.)*sin(dE3*c_PI))/(dE2*dE3*c_PI);
			a_A_avg(2) = (-4*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI)*sin(2*E3*c_PI))/(dE2*dE3*c_PI);
			a_A_avg(3) = (2*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI)*sin(2*E3*c_PI))/(dE2*dE3*pow(c_PI,2));
			a_A_avg(4) = (2*cos(E2*c_PI)*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(2*E3*c_PI))/(dE2*dE3*c_PI);
			a_A_avg(5) = (4*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI))/(dE2*dE3*c_PI);
			a_A_avg(6) = (2*cos(E2*c_PI)*sin((dE2*c_PI)/2.))/(dE2*c_PI);
			a_A_avg(7) = (-2*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;
			a_A_avg(8) = 0.0;

			double det_A_avg = a_A_avg(0)*(a_A_avg(4)*a_A_avg(8) - a_A_avg(7)*a_A_avg(5)) - a_A_avg(1)*(a_A_avg(3)*a_A_avg(8) - a_A_avg(5)*a_A_avg(6)) + a_A_avg(2)*(a_A_avg(3)*a_A_avg(7) - a_A_avg(4)*a_A_avg(6));

			a_A_inv_avg(0) = (a_A_avg(4)*a_A_avg(8) - a_A_avg(5)*a_A_avg(7))/det_A_avg;
			a_A_inv_avg(1) = (a_A_avg(2)*a_A_avg(7) - a_A_avg(1)*a_A_avg(8))/det_A_avg;
			a_A_inv_avg(2) = (a_A_avg(1)*a_A_avg(5) - a_A_avg(2)*a_A_avg(4))/det_A_avg;
			a_A_inv_avg(3) = (a_A_avg(5)*a_A_avg(6) - a_A_avg(3)*a_A_avg(8))/det_A_avg;
			a_A_inv_avg(4) = (a_A_avg(0)*a_A_avg(8) - a_A_avg(2)*a_A_avg(6))/det_A_avg;
			a_A_inv_avg(5) = (a_A_avg(2)*a_A_avg(3) - a_A_avg(0)*a_A_avg(5))/det_A_avg;
			a_A_inv_avg(6) = (a_A_avg(3)*a_A_avg(7) - a_A_avg(4)*a_A_avg(6))/det_A_avg;
			a_A_inv_avg(7) = (a_A_avg(1)*a_A_avg(6) - a_A_avg(0)*a_A_avg(7))/det_A_avg;
			a_A_inv_avg(8) = (a_A_avg(0)*a_A_avg(4) - a_A_avg(1)*a_A_avg(3))/det_A_avg;
			

			a_r2rdot_avg(0) = 1.0;
			a_detA_avg(0)	=	(4*c_PI*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;
			a_A_row_mag_avg(0) = 1.0;
			a_A_row_mag_avg(1) = c_PI;
			a_A_row_mag_avg(2) = (4*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;
		}
		if (a_r_dir_turn){
			a_r2rdot_avg(0)	*=	(exp(c*((-3*dE1)/2. + E1))*(-1 + exp(c*dE1))*Rt*(3*exp(c*dE1)*pow(R0 - Rt,2) + 3*exp(c*((3*dE1)/2. + E1))*(R0 - Rt)*Rt + 3*exp((c*dE1)/2. + c*E1)*(R0 - Rt)*Rt + exp(2*c*E1)*pow(Rt,2) + exp(2*c*(dE1 + E1))*pow(Rt,2) + exp(c*(dE1 + 2*E1))*pow(Rt,2)))/(3.*dE1);
		}
		

		a_Jacobian_ave(0) = a_r2rdot_avg(0)*a_detA_avg(0);

// FACE E=1 AVERAGED VALUES 

		E1 = (a_pt[0])*dE1;
		E2 = (a_pt[1] + 0.5)*dE2;
		E3 = (a_pt[2] + 0.5)*dE3;

		if (E2 < 0.0){
			E2 = -E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}

		if (E2 > 1.0){
			E2 = 2.0 - E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}

		if (!a_r_dir_turn){
			a_r2detA_1_avg(0)	=	(4*c_PI*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;
			a_r2detAA_1_avg(0)	=	-0.5*(cos(2*E3*c_PI)*sin(dE3*c_PI)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI)))/(dE2*dE3);
			a_r2detAA_1_avg(1)	=	(c_PI*cos(2*E3*c_PI)*sin(dE2*c_PI)*sin(dE3*c_PI)*sin(2*E2*c_PI))/(dE2*dE3);
			a_r2detAA_1_avg(2)	=	(c_PI*sin(dE3*c_PI)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI))*sin(2*E3*c_PI))/(dE2*dE3);
			a_r2detAA_1_avg(3)	=	-0.5*(sin(dE3*c_PI)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI))*sin(2*E3*c_PI))/(dE2*dE3);
			a_r2detAA_1_avg(4)	=	(c_PI*sin(dE2*c_PI)*sin(dE3*c_PI)*sin(2*E2*c_PI)*sin(2*E3*c_PI))/(dE2*dE3);
			a_r2detAA_1_avg(5)	=	-((c_PI*cos(2*E3*c_PI)*sin(dE3*c_PI)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI)))/(dE2*dE3));
			a_r2detAA_1_avg(6)	=	(c_PI*sin(dE2*c_PI)*sin(2*E2*c_PI))/dE2;
			a_r2detAA_1_avg(7)	=	(pow(c_PI,2)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI)))/(2.*dE2);
			a_r2detAA_1_avg(8)	=	0.0;
			a_r2detAn_1_avg(0)	=	-0.5*(cos(2*E3*c_PI)*sin(dE3*c_PI)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI)))/(dE2*dE3);
			a_r2detAn_1_avg(1)	=	-0.5*(sin(dE3*c_PI)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI))*sin(2*E3*c_PI))/(dE2*dE3);
			a_r2detAn_1_avg(2)	=	(c_PI*sin(dE2*c_PI)*sin(2*E2*c_PI))/dE2;

			a_n_1_avg(0)	=	(2*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI))/(dE2*dE3*pow(c_PI,2));
			a_n_1_avg(1)	=	(2*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI)*sin(2*E3*c_PI))/(dE2*dE3*pow(c_PI,2));
			a_n_1_avg(2)	=	(2*cos(E2*c_PI)*sin((dE2*c_PI)/2.))/(dE2*c_PI);

			a_A_row_mag_1_avg(0) = 1.0;
			a_A_row_mag_1_avg(1) = c_PI;
			a_A_row_mag_1_avg(2) = (4*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;

			a_A_1_avg(0) = (2*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI))/(dE2*dE3*pow(c_PI,2));
			a_A_1_avg(1) = (2*cos(E2*c_PI)*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.)*sin(dE3*c_PI))/(dE2*dE3*c_PI);
			a_A_1_avg(2) = (-4*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI)*sin(2*E3*c_PI))/(dE2*dE3*c_PI);
			a_A_1_avg(3) = (2*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI)*sin(2*E3*c_PI))/(dE2*dE3*pow(c_PI,2));
			a_A_1_avg(4) = (2*cos(E2*c_PI)*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(2*E3*c_PI))/(dE2*dE3*c_PI);
			a_A_1_avg(5) = (4*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.)*sin(dE3*c_PI)*sin(E2*c_PI))/(dE2*dE3*c_PI);
			a_A_1_avg(6) = (2*cos(E2*c_PI)*sin((dE2*c_PI)/2.))/(dE2*c_PI);
			a_A_1_avg(7) = (-2*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;
			a_A_1_avg(8) = 0.0;


			double det_A_1_avg = a_A_1_avg(0)*(a_A_1_avg(4)*a_A_1_avg(8) - a_A_1_avg(7)*a_A_1_avg(5)) - a_A_1_avg(1)*(a_A_1_avg(3)*a_A_1_avg(8) - a_A_1_avg(5)*a_A_1_avg(6)) + a_A_1_avg(2)*(a_A_1_avg(3)*a_A_1_avg(7) - a_A_1_avg(4)*a_A_1_avg(6));

			a_A_inv_1_avg(0) = (a_A_1_avg(4)*a_A_1_avg(8) - a_A_1_avg(5)*a_A_1_avg(7))/det_A_1_avg;
			a_A_inv_1_avg(1) = (a_A_1_avg(2)*a_A_1_avg(7) - a_A_1_avg(1)*a_A_1_avg(8))/det_A_1_avg;
			a_A_inv_1_avg(2) = (a_A_1_avg(1)*a_A_1_avg(5) - a_A_1_avg(2)*a_A_1_avg(4))/det_A_1_avg;
			a_A_inv_1_avg(3) = (a_A_1_avg(5)*a_A_1_avg(6) - a_A_1_avg(3)*a_A_1_avg(8))/det_A_1_avg;
			a_A_inv_1_avg(4) = (a_A_1_avg(0)*a_A_1_avg(8) - a_A_1_avg(2)*a_A_1_avg(6))/det_A_1_avg;
			a_A_inv_1_avg(5) = (a_A_1_avg(2)*a_A_1_avg(3) - a_A_1_avg(0)*a_A_1_avg(5))/det_A_1_avg;
			a_A_inv_1_avg(6) = (a_A_1_avg(3)*a_A_1_avg(7) - a_A_1_avg(4)*a_A_1_avg(6))/det_A_1_avg;
			a_A_inv_1_avg(7) = (a_A_1_avg(1)*a_A_1_avg(6) - a_A_1_avg(0)*a_A_1_avg(7))/det_A_1_avg;
			a_A_inv_1_avg(8) = (a_A_1_avg(0)*a_A_1_avg(4) - a_A_1_avg(1)*a_A_1_avg(3))/det_A_1_avg;
		}
		if (a_r_dir_turn){
			a_r2detA_1_avg(0)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAA_1_avg(0)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAA_1_avg(1)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAA_1_avg(2)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAA_1_avg(3)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAA_1_avg(4)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAA_1_avg(5)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAA_1_avg(6)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAA_1_avg(7)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAA_1_avg(8)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAn_1_avg(0)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAn_1_avg(1)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
			a_r2detAn_1_avg(2)	*=	pow(R0 + (-1 + exp(c*E1))*Rt,2);
		}

		


// FACE E=2 AVERAGED VALUES
		E1 = (a_pt[0] + 0.5)*dE1;
		E2 = (a_pt[1])*dE2;
		E3 = (a_pt[2] + 0.5)*dE3;


		if (E2 < 0.0){
			E2 = -E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}

		if (E2 > 1.0){
			E2 = 2.0 - E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}
		

		if (!a_r_dir_turn){
			a_rrdotdetA_2_avg(0)	=	(4*pow(c_PI,2)*sin(E2*c_PI));
			a_rrdotdetAA_2_avg(0)	=	(2*c_PI*pow(sin(E2*c_PI),2)*(sin((dE3 - 2*E3)*c_PI) + sin((dE3 + 2*E3)*c_PI)))/(dE3);
			a_rrdotdetAA_2_avg(1)	=	(pow(c_PI,2)*sin(2*E2*c_PI)*(sin((dE3 - 2*E3)*c_PI) + sin((dE3 + 2*E3)*c_PI)))/(dE3);
			a_rrdotdetAA_2_avg(2)	=	(-8*pow(c_PI,2)*sin(dE3*c_PI)*pow(sin(E2*c_PI),2)*sin(2*E3*c_PI))/(dE3);
			a_rrdotdetAA_2_avg(3)	=	(4*c_PI*sin(dE3*c_PI)*pow(sin(E2*c_PI),2)*sin(2*E3*c_PI))/(dE3);
			a_rrdotdetAA_2_avg(4)	=	(2*pow(c_PI,2)*sin(dE3*c_PI)*sin(2*E2*c_PI)*sin(2*E3*c_PI))/(dE3);
			a_rrdotdetAA_2_avg(5)	=	(4*pow(c_PI,2)*pow(sin(E2*c_PI),2)*(sin((dE3 - 2*E3)*c_PI) + sin((dE3 + 2*E3)*c_PI)))/(dE3);
			a_rrdotdetAA_2_avg(6)	=	(2*pow(c_PI,2)*sin(2*E2*c_PI));
			a_rrdotdetAA_2_avg(7)	=	(-4*pow(c_PI,3)*pow(sin(E2*c_PI),2));
			a_rrdotdetAA_2_avg(8)	=	0.0;
			a_rrdotd3ncn_2_avg(0)	=	(sin(2*E2*c_PI)*(sin((dE3 - 2*E3)*c_PI) + sin((dE3 + 2*E3)*c_PI)))/(dE3);
			a_rrdotd3ncn_2_avg(1)	=	(2*sin(dE3*c_PI)*sin(2*E2*c_PI)*sin(2*E3*c_PI))/(dE3);
			a_rrdotd3ncn_2_avg(2)	=	(-4*c_PI*pow(sin(E2*c_PI),2));
			a_A_row_mag_2_avg(0) = 1.0;
			a_A_row_mag_2_avg(1) = c_PI;
			a_A_row_mag_2_avg(2) = 2*c_PI*sin(E2*c_PI);

			a_A_2_avg(0) = (cos(2*E3*c_PI)*sin(dE3*c_PI)*sin(E2*c_PI))/(dE3*c_PI);
			a_A_2_avg(1) = (cos(E2*c_PI)*(sin((dE3 - 2*E3)*c_PI) + sin((dE3 + 2*E3)*c_PI)))/(2.*dE3);
			a_A_2_avg(2) = (-2*sin(dE3*c_PI)*sin(E2*c_PI)*sin(2*E3*c_PI))/dE3;
			a_A_2_avg(3) = (sin(dE3*c_PI)*sin(E2*c_PI)*sin(2*E3*c_PI))/(dE3*c_PI);
			a_A_2_avg(4) = (cos(E2*c_PI)*sin(dE3*c_PI)*sin(2*E3*c_PI))/dE3;
			a_A_2_avg(5) = (sin(E2*c_PI)*(sin((dE3 - 2*E3)*c_PI) + sin((dE3 + 2*E3)*c_PI)))/dE3;
			a_A_2_avg(6) = cos(E2*c_PI);
			a_A_2_avg(7) = -(c_PI*sin(E2*c_PI));
			a_A_2_avg(8) = 0.0;

			double det_A_2_avg = a_A_2_avg(0)*(a_A_2_avg(4)*a_A_2_avg(8) - a_A_2_avg(7)*a_A_2_avg(5)) - a_A_2_avg(1)*(a_A_2_avg(3)*a_A_2_avg(8) - a_A_2_avg(5)*a_A_2_avg(6)) + a_A_2_avg(2)*(a_A_2_avg(3)*a_A_2_avg(7) - a_A_2_avg(4)*a_A_2_avg(6));

			a_A_inv_2_avg(0) = (a_A_2_avg(4)*a_A_2_avg(8) - a_A_2_avg(5)*a_A_2_avg(7))/det_A_2_avg;
			a_A_inv_2_avg(1) = (a_A_2_avg(2)*a_A_2_avg(7) - a_A_2_avg(1)*a_A_2_avg(8))/det_A_2_avg;
			a_A_inv_2_avg(2) = (a_A_2_avg(1)*a_A_2_avg(5) - a_A_2_avg(2)*a_A_2_avg(4))/det_A_2_avg;
			a_A_inv_2_avg(3) = (a_A_2_avg(5)*a_A_2_avg(6) - a_A_2_avg(3)*a_A_2_avg(8))/det_A_2_avg;
			a_A_inv_2_avg(4) = (a_A_2_avg(0)*a_A_2_avg(8) - a_A_2_avg(2)*a_A_2_avg(6))/det_A_2_avg;
			a_A_inv_2_avg(5) = (a_A_2_avg(2)*a_A_2_avg(3) - a_A_2_avg(0)*a_A_2_avg(5))/det_A_2_avg;
			a_A_inv_2_avg(6) = (a_A_2_avg(3)*a_A_2_avg(7) - a_A_2_avg(4)*a_A_2_avg(6))/det_A_2_avg;
			a_A_inv_2_avg(7) = (a_A_2_avg(1)*a_A_2_avg(6) - a_A_2_avg(0)*a_A_2_avg(7))/det_A_2_avg;
			a_A_inv_2_avg(8) = (a_A_2_avg(0)*a_A_2_avg(4) - a_A_2_avg(1)*a_A_2_avg(3))/det_A_2_avg;
		}
		if (a_r_dir_turn){
			a_rrdotdetA_2_avg(0)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_2_avg(0)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_2_avg(1)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_2_avg(2)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_2_avg(3)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_2_avg(4)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_2_avg(5)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_2_avg(6)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_2_avg(7)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_2_avg(8)	*=	0.0;
			a_rrdotd3ncn_2_avg(0)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotd3ncn_2_avg(1)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotd3ncn_2_avg(2)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
		}

		


// FACE E=3 AVERAGED VALUES

		E1 = (a_pt[0] + 0.5)*dE1;
		E2 = (a_pt[1] + 0.5)*dE2;
		E3 = (a_pt[2])*dE3;

		
		if (E2 < 0.0){
			E2 = -E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}

		if (E2 > 1.0){
			E2 = 2.0 - E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}
		
		
		if (!a_r_dir_turn){
			a_rrdotdetA_3_avg(0)	=	(8*c_PI*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/(dE2);
			a_rrdotdetAA_3_avg(0)	=	-((cos(2*E3*c_PI)*c_PI*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI)))/(dE2));
			a_rrdotdetAA_3_avg(1)	=	(2*pow(c_PI,2)*cos(2*E3*c_PI)*sin(dE2*c_PI)*sin(2*E2*c_PI))/(dE2);
			a_rrdotdetAA_3_avg(2)	=	(2*pow(c_PI,2)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI))*sin(2*E3*c_PI))/(dE2);
			a_rrdotdetAA_3_avg(3)	=	-((c_PI*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI))*sin(2*E3*c_PI))/(dE2));
			a_rrdotdetAA_3_avg(4)	=	(2*pow(c_PI,2)*sin(dE2*c_PI)*sin(2*E2*c_PI)*sin(2*E3*c_PI))/(dE2);
			a_rrdotdetAA_3_avg(5)	=	(-2*pow(c_PI,2)*cos(2*E3*c_PI)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI)))/(dE2);
			a_rrdotdetAA_3_avg(6)	=	(2*c_PI*sin(dE2*c_PI)*sin(2*E2*c_PI))/(dE2);
			a_rrdotdetAA_3_avg(7)	=	(pow(c_PI,2)*(-2*dE2*c_PI + sin((dE2 - 2*E2)*c_PI) + sin((dE2 + 2*E2)*c_PI)))/(dE2);
			a_rrdotdetAA_3_avg(8)	=	0.0;
			a_rrdotncd2n_3_avg(0)	=	(-2*c_PI*sin(2*E3*c_PI));
			a_rrdotncd2n_3_avg(1)	=	(2*c_PI*cos(2*E3*c_PI));
			a_rrdotncd2n_3_avg(2)	=	0.0;
			a_A_row_mag_3_avg(0) = 1.0;
			a_A_row_mag_3_avg(1) = c_PI;
			a_A_row_mag_3_avg(2) = (4*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;

			a_A_3_avg(0) = (2*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/(dE2*c_PI);
			a_A_3_avg(1) = (2*cos(E2*c_PI)*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.))/dE2;
			a_A_3_avg(2) = (-4*sin((dE2*c_PI)/2.)*sin(E2*c_PI)*sin(2*E3*c_PI))/dE2;
			a_A_3_avg(3) = (2*sin((dE2*c_PI)/2.)*sin(E2*c_PI)*sin(2*E3*c_PI))/(dE2*c_PI);
			a_A_3_avg(4) = (2*cos(E2*c_PI)*sin((dE2*c_PI)/2.)*sin(2*E3*c_PI))/dE2;
			a_A_3_avg(5) = (4*cos(2*E3*c_PI)*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;
			a_A_3_avg(6) = (2*cos(E2*c_PI)*sin((dE2*c_PI)/2.))/(dE2*c_PI);
			a_A_3_avg(7) = (-2*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;
			a_A_3_avg(8) = 0.0;

			double det_A_3_avg = a_A_3_avg(0)*(a_A_3_avg(4)*a_A_3_avg(8) - a_A_3_avg(7)*a_A_3_avg(5)) - a_A_3_avg(1)*(a_A_3_avg(3)*a_A_3_avg(8) - a_A_3_avg(5)*a_A_3_avg(6)) + a_A_3_avg(2)*(a_A_3_avg(3)*a_A_3_avg(7) - a_A_3_avg(4)*a_A_3_avg(6));

			a_A_inv_3_avg(0) = (a_A_3_avg(4)*a_A_3_avg(8) - a_A_3_avg(5)*a_A_3_avg(7))/det_A_3_avg;
			a_A_inv_3_avg(1) = (a_A_3_avg(2)*a_A_3_avg(7) - a_A_3_avg(1)*a_A_3_avg(8))/det_A_3_avg;
			a_A_inv_3_avg(2) = (a_A_3_avg(1)*a_A_3_avg(5) - a_A_3_avg(2)*a_A_3_avg(4))/det_A_3_avg;
			a_A_inv_3_avg(3) = (a_A_3_avg(5)*a_A_3_avg(6) - a_A_3_avg(3)*a_A_3_avg(8))/det_A_3_avg;
			a_A_inv_3_avg(4) = (a_A_3_avg(0)*a_A_3_avg(8) - a_A_3_avg(2)*a_A_3_avg(6))/det_A_3_avg;
			a_A_inv_3_avg(5) = (a_A_3_avg(2)*a_A_3_avg(3) - a_A_3_avg(0)*a_A_3_avg(5))/det_A_3_avg;
			a_A_inv_3_avg(6) = (a_A_3_avg(3)*a_A_3_avg(7) - a_A_3_avg(4)*a_A_3_avg(6))/det_A_3_avg;
			a_A_inv_3_avg(7) = (a_A_3_avg(1)*a_A_3_avg(6) - a_A_3_avg(0)*a_A_3_avg(7))/det_A_3_avg;
			a_A_inv_3_avg(8) = (a_A_3_avg(0)*a_A_3_avg(4) - a_A_3_avg(1)*a_A_3_avg(3))/det_A_3_avg;
		}
		if (a_r_dir_turn){
			a_rrdotdetA_3_avg(0)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_3_avg(0)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_3_avg(1)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_3_avg(2)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_3_avg(3)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_3_avg(4)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_3_avg(5)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_3_avg(6)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_3_avg(7)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotdetAA_3_avg(8)	*=	0.0;
			a_rrdotncd2n_3_avg(0)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotncd2n_3_avg(1)	*=	(exp(c*E1)*Rt*(R0 - Rt + exp(c*E1)*Rt*cosh((c*dE1)/2.))*sinh((c*dE1)/2.))/dE1;
			a_rrdotncd2n_3_avg(2)	*=	0.0;
		}

		

	}
	PROTO_KERNEL_END(Spherical_map_calcF, Spherical_map_calc)



	void Spherical_map_calc_func(BoxData<double,1>& a_Jacobian_ave,
								 BoxData<double,DIM*DIM>& a_A_avg,
								 BoxData<double,DIM*DIM>& a_A_inv_avg,
	                             BoxData<double,DIM*DIM>& a_A_1_avg,
	                             BoxData<double,DIM*DIM>& a_A_2_avg,
	                             BoxData<double,DIM*DIM>& a_A_3_avg,
								 BoxData<double,DIM*DIM>& a_A_inv_1_avg,
	                             BoxData<double,DIM*DIM>& a_A_inv_2_avg,
	                             BoxData<double,DIM*DIM>& a_A_inv_3_avg,
	                             BoxData<double,DIM*DIM>& a_detAA_avg,
	                             BoxData<double,DIM*DIM>& a_detAA_inv_avg,
	                             BoxData<double,1>& a_r2rdot_avg,
	                             BoxData<double,1>& a_detA_avg,
	                             BoxData<double,DIM>& a_A_row_mag_avg,
	                             BoxData<double,1>& a_r2detA_1_avg,
	                             BoxData<double,DIM*DIM>& a_r2detAA_1_avg,
	                             BoxData<double,DIM>& a_r2detAn_1_avg,
	                             BoxData<double,DIM>& a_n_1_avg,
								 BoxData<double,DIM>& a_A_row_mag_1_avg,
	                             BoxData<double,1>& a_rrdotdetA_2_avg,
	                             BoxData<double,DIM*DIM>& a_rrdotdetAA_2_avg,
	                             BoxData<double,DIM>& a_rrdotd3ncn_2_avg,
								 BoxData<double,DIM>& a_A_row_mag_2_avg,
	                             BoxData<double,1>& a_rrdotdetA_3_avg,
	                             BoxData<double,DIM*DIM>& a_rrdotdetAA_3_avg,
	                             BoxData<double,DIM>& a_rrdotncd2n_3_avg,
								 BoxData<double,DIM>& a_A_row_mag_3_avg,
	                             const double a_dx,
	                             const double a_dy,
	                             const double a_dz,
								 bool a_exchanged_yet,
								 bool a_r_dir_turn)
	{
		forallInPlace_p(Spherical_map_calc, a_Jacobian_ave, a_A_avg, a_A_inv_avg, a_A_1_avg, a_A_2_avg, a_A_3_avg, a_A_inv_1_avg, a_A_inv_2_avg, a_A_inv_3_avg, a_detAA_avg, a_detAA_inv_avg, a_r2rdot_avg, a_detA_avg, a_A_row_mag_avg, a_r2detA_1_avg, a_r2detAA_1_avg, a_r2detAn_1_avg, a_n_1_avg, a_A_row_mag_1_avg, a_rrdotdetA_2_avg, a_rrdotdetAA_2_avg, a_rrdotd3ncn_2_avg, a_A_row_mag_2_avg, a_rrdotdetA_3_avg, a_rrdotdetAA_3_avg, a_rrdotncd2n_3_avg, a_A_row_mag_3_avg, a_dx, a_dy, a_dz,a_exchanged_yet,a_r_dir_turn);
	}

	PROTO_KERNEL_START
	void Jacobian_ave_sph_calcF( const Point& a_pt,
	                          Var<double,1>& a_Jacobian_ave,
	                          const double a_dx,
	                          const double a_dy,
	                          const double a_dz)
	{

		double R0 = inputs.r_in*c_AU;
		double R1 = inputs.r_out*c_AU;
		double c = inputs.C_rad;
		double Rt = (R1 - R0)/(exp(c) - 1.0);

		double dE1 = a_dx;
		double dE2 = a_dy;
		double dE3 = a_dz;

		double E1, E2, E3;
/////  CELL AVERAGED VALUES

		//For the special mapping that preserves radial flow, Phil has suggested to use E2 as theta and E3 as Phi.

		E1 = (a_pt[0] + 0.5)*dE1;
		E2 = (a_pt[1] + 0.5)*dE2;
		E3 = (a_pt[2] + 0.5)*dE3;

		
		if (E2 < 0.0){
			E2 = -E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}

		if (E2 > 1.0){
			E2 = 2.0 - E2;
			if (E3 < 0.5) {
				E3 += 0.5;
			} else {
				E3 -= 0.5;
			}
		}
		
		double a_r2rdot_avg	=	(exp(c*((-3*dE1)/2. + E1))*(-1 + exp(c*dE1))*Rt*(3*exp(c*dE1)*pow(R0 - Rt,2) + 3*exp(c*((3*dE1)/2. + E1))*(R0 - Rt)*Rt + 3*exp((c*dE1)/2. + c*E1)*(R0 - Rt)*Rt + exp(2*c*E1)*pow(Rt,2) + exp(2*c*(dE1 + E1))*pow(Rt,2) + exp(c*(dE1 + 2*E1))*pow(Rt,2)))/(3.*dE1);
		double a_detA_avg	=	(4*c_PI*sin((dE2*c_PI)/2.)*sin(E2*c_PI))/dE2;

		a_Jacobian_ave(0) = a_r2rdot_avg*a_detA_avg;
	}
	PROTO_KERNEL_END(Jacobian_ave_sph_calcF, Jacobian_ave_sph_calc)


	void Jacobian_ave_sph_calc_func(BoxData<double,1>& a_Jacobian_ave,
	                             const double a_dx,
	                             const double a_dy,
	                             const double a_dz)
	{
		forallInPlace_p(Jacobian_ave_sph_calc, a_Jacobian_ave, a_dx, a_dy, a_dz);
	}


	PROTO_KERNEL_START
	void JU_to_U_Sph_ave_calcF(const Point& a_pt,
							   State& a_U_Sph_ave,
	                          const Var<double,NUMCOMPS>& a_JU_ave,
	                          const Var<double,DIM*DIM>& a_detAA_inv_avg,
	                          const Var<double,DIM*DIM>& a_A_inv_avg,
	                          const Var<double,1>& a_r2rdot_avg,
	                          const Var<double,1>& a_detA_avg,
	                          const Var<double,DIM>& a_A_row_mag_avg,
							  bool a_normalized)
	{
		double r2rdotrho_ave = a_JU_ave(0)/a_detA_avg(0);
		a_U_Sph_ave(0) = r2rdotrho_ave/a_r2rdot_avg(0);

		double r2rdotE_ave = a_JU_ave(4)/a_detA_avg(0);
		a_U_Sph_ave(4) = r2rdotE_ave/a_r2rdot_avg(0);

		double r2rdotrhou_ave = a_detAA_inv_avg(0)*a_JU_ave(1) + a_detAA_inv_avg(1)*a_JU_ave(2) + a_detAA_inv_avg(2)*a_JU_ave(3);
		double r2rdotrhov_ave = a_detAA_inv_avg(3)*a_JU_ave(1) + a_detAA_inv_avg(4)*a_JU_ave(2) + a_detAA_inv_avg(5)*a_JU_ave(3);
		double r2rdotrhow_ave = a_detAA_inv_avg(6)*a_JU_ave(1) + a_detAA_inv_avg(7)*a_JU_ave(2) + a_detAA_inv_avg(8)*a_JU_ave(3);

		// double r2rdotrhou_ave = a_A_inv_avg(0)*a_JU_ave(1) + a_A_inv_avg(1)*a_JU_ave(2) + a_A_inv_avg(2)*a_JU_ave(3);
		// double r2rdotrhov_ave = a_A_inv_avg(3)*a_JU_ave(1) + a_A_inv_avg(4)*a_JU_ave(2) + a_A_inv_avg(5)*a_JU_ave(3);
		// double r2rdotrhow_ave = a_A_inv_avg(6)*a_JU_ave(1) + a_A_inv_avg(7)*a_JU_ave(2) + a_A_inv_avg(8)*a_JU_ave(3);

		a_U_Sph_ave(1)  = r2rdotrhou_ave/a_r2rdot_avg(0);
		a_U_Sph_ave(2)  = r2rdotrhov_ave/a_r2rdot_avg(0);
		a_U_Sph_ave(3)  = r2rdotrhow_ave/a_r2rdot_avg(0);

		// a_U_Sph_ave(1)  = r2rdotrhou_ave/a_r2rdot_avg(0)/a_detA_avg(0);
		// a_U_Sph_ave(2)  = r2rdotrhov_ave/a_r2rdot_avg(0)/a_detA_avg(0);
		// a_U_Sph_ave(3)  = r2rdotrhow_ave/a_r2rdot_avg(0)/a_detA_avg(0);

		double r2rdotBx_ave = a_detAA_inv_avg(0)*a_JU_ave(5) + a_detAA_inv_avg(1)*a_JU_ave(6) + a_detAA_inv_avg(2)*a_JU_ave(7);
		double r2rdotBy_ave = a_detAA_inv_avg(3)*a_JU_ave(5) + a_detAA_inv_avg(4)*a_JU_ave(6) + a_detAA_inv_avg(5)*a_JU_ave(7);
		double r2rdotBz_ave = a_detAA_inv_avg(6)*a_JU_ave(5) + a_detAA_inv_avg(7)*a_JU_ave(6) + a_detAA_inv_avg(8)*a_JU_ave(7);
		a_U_Sph_ave(5)  = r2rdotBx_ave/a_r2rdot_avg(0);
		a_U_Sph_ave(6)  = r2rdotBy_ave/a_r2rdot_avg(0);
		a_U_Sph_ave(7)  = r2rdotBz_ave/a_r2rdot_avg(0);

		if (a_normalized){
			a_U_Sph_ave(1)*=a_A_row_mag_avg(0);
			a_U_Sph_ave(2)*=a_A_row_mag_avg(1);
			a_U_Sph_ave(3)*=a_A_row_mag_avg(2);
			a_U_Sph_ave(5)*=a_A_row_mag_avg(0);
			a_U_Sph_ave(6)*=a_A_row_mag_avg(1);
			a_U_Sph_ave(7)*=a_A_row_mag_avg(2);
		}
	}
	PROTO_KERNEL_END(JU_to_U_Sph_ave_calcF, JU_to_U_Sph_ave_calc)

	void JU_to_U_Sph_ave_calc_func(BoxData<double,NUMCOMPS>& a_U_Sph_ave,
	                  const BoxData<double,NUMCOMPS>& a_JU_ave,
	                  BoxData<double,DIM*DIM>& a_detAA_inv_avg,
	                  BoxData<double,DIM*DIM>& a_A_inv_avg,
	                  BoxData<double,1>& a_r2rdot_avg,
	                  BoxData<double,1>& a_detA_avg,
	                  BoxData<double,DIM>& a_A_row_mag_avg,
					  bool a_normalized,
					  int a_order)
	{
		if (a_order == 2){
			forallInPlace_p(JU_to_U_Sph_ave_calc, a_U_Sph_ave, a_JU_ave,a_detAA_inv_avg, a_A_inv_avg,a_r2rdot_avg,a_detA_avg, a_A_row_mag_avg, a_normalized);
		}
		if (a_order == 4){
			#define CRHO 0
			#define CVELSTART 1
			#define CBSTART 5
			#define CPRES 4
			#define CENG 4
			a_U_Sph_ave.setVal(0.0);
			Box dbx0 = a_JU_ave.box();
			cout << dbx0 << endl;
			BoxData<double,DIM,MEM,DIM> Ainv4(dbx0);
			MHD_Mapping::Nineto33(Ainv4, a_detAA_inv_avg);
			auto JU4 = slice<double,NUMCOMPS,DIM,MEM>(a_JU_ave,CVELSTART);
			auto JB4 = slice<double,NUMCOMPS,DIM,MEM>(a_JU_ave,CBSTART);  
			auto JU2 = slice<double,NUMCOMPS,DIM,MEM>(a_JU_ave,CVELSTART);
			auto JB2 = slice<double,NUMCOMPS,DIM,MEM>(a_JU_ave,CBSTART);
			BoxData<double,DIM,MEM> w4(dbx0), b4(dbx0), U4(dbx0), B4(dbx0);
			BoxData<double,DIM,MEM> U4_temp = Operator::_matrixProductAB2(Ainv4,JU4);
			BoxData<double,DIM,MEM> B4_temp = Operator::_matrixProductAB2(Ainv4,JB4);
			if (a_normalized){
				w4 = Operator::_cellTensorQuotient(U4_temp,a_r2rdot_avg,U4_temp,a_r2rdot_avg);
				b4 = Operator::_cellTensorQuotient(B4_temp,a_r2rdot_avg,B4_temp,a_r2rdot_avg);
				U4 = Operator::cellProduct(w4,a_A_row_mag_avg);
				B4 = Operator::cellProduct(b4,a_A_row_mag_avg);
			} else {
				U4 = Operator::_cellTensorQuotient(U4_temp,a_r2rdot_avg,U4_temp,a_r2rdot_avg);
				B4 = Operator::_cellTensorQuotient(B4_temp,a_r2rdot_avg,B4_temp,a_r2rdot_avg);
			}
			
			
			BoxData<double,1,MEM> Jrho4 = slice(a_JU_ave,CRHO);
			auto rho4_temp = Operator::_cellQuotient(Jrho4,a_detA_avg,Jrho4,a_detA_avg);
			auto rho4 = Operator::_cellQuotient(rho4_temp,a_r2rdot_avg,rho4_temp,a_r2rdot_avg);

			BoxData<double,1,MEM> JE4 = slice(a_JU_ave,CENG);
			auto E4_temp = Operator::cellQuotient(JE4,a_detA_avg);
			auto E4 = Operator::cellQuotient(E4_temp,a_r2rdot_avg);

			a_U_Sph_ave = forall<double,NUMCOMPS,MEM,1>
				([ ] PROTO_LAMBDA(
								Var<double,NUMCOMPS,MEM,1>& a_retval,
								Var<double,1,MEM>& a_rho4,
								Var<double,DIM,MEM>& a_U4,
								Var<double,1,MEM>& a_E4,
								Var<double,DIM,MEM>& a_B4)
				{
				a_retval(0) = a_rho4(0);
				a_retval(4) = a_E4(0);
				for (int dir = 0; dir < DIM; dir++)
					{
					a_retval(1+dir) = a_U4(dir);
					a_retval(5+dir) = a_B4(dir);
					}
				},
				rho4,U4,E4,B4);
		}
	}

	void JU_to_W_Sph_ave_calc_func(BoxData<double,NUMCOMPS>& a_W_Sph_ave,
	                  const BoxData<double,NUMCOMPS>& a_JU_ave,
	                  BoxData<double,DIM*DIM>& a_detAA_inv_avg,
	                  BoxData<double,DIM*DIM>& a_A_inv_avg,
	                  BoxData<double,1>& a_r2rdot_avg,
	                  BoxData<double,1>& a_detA_avg,
					  BoxData<double,DIM>& a_A_row_mag_avg,
	                  const double a_gamma,
					  bool a_normalized)
	{
		double gamma = a_gamma;
		Vector a_U_Sph_ave(a_W_Sph_ave.box());
		MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_ave, a_JU_ave,a_detAA_inv_avg,a_A_inv_avg,a_r2rdot_avg,a_detA_avg, a_A_row_mag_avg, a_normalized, 2);
		MHDOp::consToPrimcalc(a_W_Sph_ave,a_U_Sph_ave,gamma);
	}



	PROTO_KERNEL_START
	void U_Sph_ave_to_JU_calcF(State& a_JU_ave,
	                          const Var<double,NUMCOMPS>& a_U_Sph_ave,
	                          const Var<double,DIM*DIM>& a_detAA_avg,
	                          const Var<double,DIM*DIM>& a_A_avg,
	                          const Var<double,DIM*DIM>& a_detAA_inv_avg,
	                          const Var<double,1>& a_r2rdot_avg,
	                          const Var<double,1>& a_detA_avg,
							  const Var<double,DIM>& a_A_row_mag_avg,
					  		  bool a_normalized)
	{

		double r2rdotrho_ave = a_U_Sph_ave(0)*a_r2rdot_avg(0);
		a_JU_ave(0) = r2rdotrho_ave * a_detA_avg(0);

		double r2rdotE_ave = a_U_Sph_ave(4)*a_r2rdot_avg(0);
		a_JU_ave(4) = r2rdotE_ave * a_detA_avg(0);

		double r2rdotrhou_ave = a_U_Sph_ave(1)*a_r2rdot_avg(0);
		double r2rdotrhov_ave = a_U_Sph_ave(2)*a_r2rdot_avg(0);
		double r2rdotrhow_ave = a_U_Sph_ave(3)*a_r2rdot_avg(0);

		double r2rdotBx_ave = a_U_Sph_ave(5)*a_r2rdot_avg(0);
		double r2rdotBy_ave = a_U_Sph_ave(6)*a_r2rdot_avg(0);
		double r2rdotBz_ave = a_U_Sph_ave(7)*a_r2rdot_avg(0);

		if (a_normalized){
			r2rdotrhou_ave/=a_A_row_mag_avg(0);
			r2rdotrhov_ave/=a_A_row_mag_avg(1);
			r2rdotrhow_ave/=a_A_row_mag_avg(2);

			r2rdotBx_ave/=a_A_row_mag_avg(0);
			r2rdotBy_ave/=a_A_row_mag_avg(1);
			r2rdotBz_ave/=a_A_row_mag_avg(2);
		}

		a_JU_ave(1) = a_detAA_avg(0)*r2rdotrhou_ave + a_detAA_avg(1)*r2rdotrhov_ave + a_detAA_avg(2)*r2rdotrhow_ave;
		a_JU_ave(2) = a_detAA_avg(3)*r2rdotrhou_ave + a_detAA_avg(4)*r2rdotrhov_ave + a_detAA_avg(5)*r2rdotrhow_ave;
		a_JU_ave(3) = a_detAA_avg(6)*r2rdotrhou_ave + a_detAA_avg(7)*r2rdotrhov_ave + a_detAA_avg(8)*r2rdotrhow_ave;

		// a_JU_ave(1) = a_detA_avg(0)*(a_A_avg(0)*r2rdotrhou_ave + a_A_avg(1)*r2rdotrhov_ave + a_A_avg(2)*r2rdotrhow_ave);
		// a_JU_ave(2) = a_detA_avg(0)*(a_A_avg(3)*r2rdotrhou_ave + a_A_avg(4)*r2rdotrhov_ave + a_A_avg(5)*r2rdotrhow_ave);
		// a_JU_ave(3) = a_detA_avg(0)*(a_A_avg(6)*r2rdotrhou_ave + a_A_avg(7)*r2rdotrhov_ave + a_A_avg(8)*r2rdotrhow_ave);

		a_JU_ave(5) = a_detAA_avg(0)*r2rdotBx_ave + a_detAA_avg(1)*r2rdotBy_ave + a_detAA_avg(2)*r2rdotBz_ave;
		a_JU_ave(6) = a_detAA_avg(3)*r2rdotBx_ave + a_detAA_avg(4)*r2rdotBy_ave + a_detAA_avg(5)*r2rdotBz_ave;
		a_JU_ave(7) = a_detAA_avg(6)*r2rdotBx_ave + a_detAA_avg(7)*r2rdotBy_ave + a_detAA_avg(8)*r2rdotBz_ave;
	}
	PROTO_KERNEL_END(U_Sph_ave_to_JU_calcF, U_Sph_ave_to_JU_calc)

	//With Phil's operators
	void U_Sph_ave_to_JU_calc_func(BoxData<double,NUMCOMPS>& a_JU_ave,
	                  const BoxData<double,NUMCOMPS>& a_U_Sph4_ave,
	                  BoxData<double,DIM*DIM>& a_detAA_avg,
	                  BoxData<double,DIM*DIM>& a_A_avg,
	                  BoxData<double,DIM*DIM>& a_detAA_inv_avg,
	                  BoxData<double,1>& a_r2rdot_avg,
	                  BoxData<double,1>& a_detA_avg,
					  BoxData<double,DIM>& a_A_row_mag_avg,
					  bool a_normalized,
					  int order)
	{
		if (order == 2){
			a_JU_ave = forall<double,NUMCOMPS>(U_Sph_ave_to_JU_calc,a_U_Sph4_ave,a_detAA_avg, a_A_avg, a_detAA_inv_avg,a_r2rdot_avg,a_detA_avg,a_A_row_mag_avg,a_normalized);
		}

		if (order == 4){
			#define CRHO 0
			#define CVELSTART 1
			#define CBSTART 5
			#define CPRES 4
			#define CENG 4
			a_JU_ave.setVal(0.0);
			Box dbx0 = a_U_Sph4_ave.box();
			cout << dbx0 << endl;
			BoxData<double,DIM,MEM,DIM> A4(dbx0);
			MHD_Mapping::Nineto33(A4, a_detAA_avg);
			auto w4 = slice<double,NUMCOMPS,DIM,MEM>(a_U_Sph4_ave,CVELSTART);
			auto b4 = slice<double,NUMCOMPS,DIM,MEM>(a_U_Sph4_ave,CBSTART);  
			auto w2 = slice<double,NUMCOMPS,DIM,MEM>(a_U_Sph4_ave,CVELSTART);
			auto b2 = slice<double,NUMCOMPS,DIM,MEM>(a_U_Sph4_ave,CBSTART);
			BoxData<double,DIM,MEM> w4_temp(dbx0), b4_temp(dbx0), w2_temp(dbx0), b2_temp(dbx0), U4_temp(dbx0), B4_temp(dbx0), U2_temp(dbx0), B2_temp(dbx0);
			
			if (a_normalized){
				w4_temp = Operator::cellQuotient(w4,a_A_row_mag_avg);
				b4_temp = Operator::cellQuotient(b4,a_A_row_mag_avg);
				w2_temp = Operator::cellQuotient(w2,a_A_row_mag_avg);
				b2_temp = Operator::cellQuotient(b2,a_A_row_mag_avg);
				U4_temp = Operator::_cellTensorProduct(w4_temp,a_r2rdot_avg,w2_temp,a_r2rdot_avg);
				B4_temp = Operator::_cellTensorProduct(b4_temp,a_r2rdot_avg,b2_temp,a_r2rdot_avg);
			} else {
				U4_temp = Operator::_cellTensorProduct(w4,a_r2rdot_avg,w2,a_r2rdot_avg);
				B4_temp = Operator::_cellTensorProduct(b4,a_r2rdot_avg,b2,a_r2rdot_avg);
			}
			BoxData<double,DIM,MEM> JU4 = Operator::_matrixProductAB2(A4,U4_temp);
			BoxData<double,DIM,MEM> JB4 = Operator::_matrixProductAB2(A4,B4_temp);
			
			BoxData<double,1,MEM> rho4 = slice(a_U_Sph4_ave,CRHO);
			BoxData<double,1,MEM> rho2 = slice(a_U_Sph4_ave,CRHO);
			auto Jrho4_temp = Operator::_cellProduct(rho4,a_r2rdot_avg,rho2,a_r2rdot_avg);
			auto Jrho4 = Operator::_cellProduct(Jrho4_temp,a_detA_avg,Jrho4_temp,a_detA_avg);

			BoxData<double,1,MEM> E4 = slice(a_U_Sph4_ave,CENG); 
			BoxData<double,1,MEM> E2 = slice(a_U_Sph4_ave,CENG);
			BoxData<double,1,MEM> JE4_temp = Operator::_cellProduct(E4,a_r2rdot_avg,E2,a_r2rdot_avg);
			BoxData<double,1,MEM> JE4 = Operator::_cellProduct(JE4_temp,a_detA_avg,JE4_temp,a_detA_avg);

			a_JU_ave = forall<double,NUMCOMPS,MEM,1>
				([ ] PROTO_LAMBDA(
								Var<double,NUMCOMPS,MEM,1>& a_retval,
								Var<double,1,MEM>& a_Jrho4,
								Var<double,DIM,MEM>& a_JU4,
								Var<double,1,MEM>& a_JE4,
								Var<double,DIM,MEM>& a_JB4)
				{
				a_retval(0) = a_Jrho4(0);
				a_retval(4) = a_JE4(0);
				for (int dir = 0; dir < DIM; dir++)
					{
					a_retval(1+dir) = a_JU4(dir);
					a_retval(5+dir) = a_JB4(dir);
					}
				},
				Jrho4,JU4,JE4,JB4);
		}
		
	}

	PROTO_KERNEL_START
	void JU_to_U_ave_calcF(State& a_U_Sph_ave,
	                          const Var<double,NUMCOMPS>& a_JU_ave,
	                          const Var<double,1>& a_r2rdot_avg,
	                          const Var<double,1>& a_detA_avg)
	{
		double r2rdotrho_ave = a_JU_ave(0)/a_detA_avg(0);
		a_U_Sph_ave(0) = r2rdotrho_ave/a_r2rdot_avg(0);

		double r2rdotE_ave = a_JU_ave(4)/a_detA_avg(0);
		a_U_Sph_ave(4) = r2rdotE_ave/a_r2rdot_avg(0);

		double r2rdotrhou_ave = a_JU_ave(1)/a_detA_avg(0);
		double r2rdotrhov_ave = a_JU_ave(2)/a_detA_avg(0);
		double r2rdotrhow_ave = a_JU_ave(3)/a_detA_avg(0);
		a_U_Sph_ave(1)  = r2rdotrhou_ave/a_r2rdot_avg(0);
		a_U_Sph_ave(2)  = r2rdotrhov_ave/a_r2rdot_avg(0);
		a_U_Sph_ave(3)  = r2rdotrhow_ave/a_r2rdot_avg(0);

		double r2rdotBx_ave = a_JU_ave(5)/a_detA_avg(0);
		double r2rdotBy_ave = a_JU_ave(6)/a_detA_avg(0);
		double r2rdotBz_ave = a_JU_ave(7)/a_detA_avg(0);
		a_U_Sph_ave(5)  = r2rdotBx_ave/a_r2rdot_avg(0);
		a_U_Sph_ave(6)  = r2rdotBy_ave/a_r2rdot_avg(0);
		a_U_Sph_ave(7)  = r2rdotBz_ave/a_r2rdot_avg(0);
	}
	PROTO_KERNEL_END(JU_to_U_ave_calcF, JU_to_U_ave_calc)

	void JU_to_U_ave_calc_func(BoxData<double,NUMCOMPS>& a_U_ave,
	                  const BoxData<double,NUMCOMPS>& a_JU_ave,
	                  BoxData<double,1>& a_r2rdot_avg,
	                  BoxData<double,1>& a_detA_avg)
	{

		a_U_ave = forall<double,NUMCOMPS>(JU_to_U_ave_calc,a_JU_ave,a_r2rdot_avg,a_detA_avg);
	}

	PROTO_KERNEL_START
	void W_Sph_to_W_Cart_calcF( const Point& a_pt,
	                          Var<double,NUMCOMPS>& W_cart,
	                          Var<double,NUMCOMPS>& W,
							  const Var<double,DIM*DIM>& a_A_1_avg,
							  const Var<double,DIM*DIM>& a_A_2_avg,
							  const Var<double,DIM*DIM>& a_A_3_avg,
							  int a_d)
	{

		if (a_d == 0)
		{
			W_cart(0) = W(0);
			W_cart(4) = W(4);
			W_cart(1) = a_A_1_avg(0)*W(1) + a_A_1_avg(1)*W(2) + a_A_1_avg(2)*W(3);
			W_cart(2) = a_A_1_avg(3)*W(1) + a_A_1_avg(4)*W(2) + a_A_1_avg(5)*W(3);
			W_cart(3) = a_A_1_avg(6)*W(1) + a_A_1_avg(7)*W(2) + a_A_1_avg(8)*W(3);
			W_cart(5) = a_A_1_avg(0)*W(5) + a_A_1_avg(1)*W(6) + a_A_1_avg(2)*W(7);
			W_cart(6) = a_A_1_avg(3)*W(5) + a_A_1_avg(4)*W(6) + a_A_1_avg(5)*W(7);
			W_cart(7) = a_A_1_avg(6)*W(5) + a_A_1_avg(7)*W(6) + a_A_1_avg(8)*W(7);
		}

		if (a_d == 1)
		{
			W_cart(0) = W(0);
			W_cart(4) = W(4);
			W_cart(1) = a_A_2_avg(0)*W(1) + a_A_2_avg(1)*W(2) + a_A_2_avg(2)*W(3);
			W_cart(2) = a_A_2_avg(3)*W(1) + a_A_2_avg(4)*W(2) + a_A_2_avg(5)*W(3);
			W_cart(3) = a_A_2_avg(6)*W(1) + a_A_2_avg(7)*W(2) + a_A_2_avg(8)*W(3);
			W_cart(5) = a_A_2_avg(0)*W(5) + a_A_2_avg(1)*W(6) + a_A_2_avg(2)*W(7);
			W_cart(6) = a_A_2_avg(3)*W(5) + a_A_2_avg(4)*W(6) + a_A_2_avg(5)*W(7);
			W_cart(7) = a_A_2_avg(6)*W(5) + a_A_2_avg(7)*W(6) + a_A_2_avg(8)*W(7);
		}

		if (a_d == 2)
		{
			W_cart(0) = W(0);
			W_cart(4) = W(4);
			W_cart(1) = a_A_3_avg(0)*W(1) + a_A_3_avg(1)*W(2) + a_A_3_avg(2)*W(3);
			W_cart(2) = a_A_3_avg(3)*W(1) + a_A_3_avg(4)*W(2) + a_A_3_avg(5)*W(3);
			W_cart(3) = a_A_3_avg(6)*W(1) + a_A_3_avg(7)*W(2) + a_A_3_avg(8)*W(3);
			W_cart(5) = a_A_3_avg(0)*W(5) + a_A_3_avg(1)*W(6) + a_A_3_avg(2)*W(7);
			W_cart(6) = a_A_3_avg(3)*W(5) + a_A_3_avg(4)*W(6) + a_A_3_avg(5)*W(7);
			W_cart(7) = a_A_3_avg(6)*W(5) + a_A_3_avg(7)*W(6) + a_A_3_avg(8)*W(7);
		}
		
	}
	PROTO_KERNEL_END(W_Sph_to_W_Cart_calcF, W_Sph_to_W_Cart_calc)

	void W_Sph_to_W_Cart(BoxData<double,NUMCOMPS>& W_cart,
	                    const BoxData<double,NUMCOMPS>& W,
						BoxData<double,DIM*DIM>& a_A_1_avg,
						BoxData<double,DIM*DIM>& a_A_2_avg,
						BoxData<double,DIM*DIM>& a_A_3_avg,
	                    int a_d)
	{
		forallInPlace_p(W_Sph_to_W_Cart_calc, W_cart, W, a_A_1_avg, a_A_2_avg, a_A_3_avg, a_d);
	}	


	PROTO_KERNEL_START
	void W_Sph_to_W_normalized_sph_calcF( const Point& a_pt,
										Var<double,NUMCOMPS>& W_normalized_sph,
										Var<double,NUMCOMPS>& W_Sph,
										const Var<double,DIM>& a_A_row_mag_1_avg,
										const Var<double,DIM>& a_A_row_mag_2_avg,
										const Var<double,DIM>& a_A_row_mag_3_avg,
										int a_d)
	{
		double a,b,c;

		if (a_d == 0)
		{
			a = a_A_row_mag_1_avg(0);
			b = a_A_row_mag_1_avg(1);
			c = a_A_row_mag_1_avg(2);
		}
		if (a_d == 1)
		{
			a = a_A_row_mag_2_avg(0);
			b = a_A_row_mag_2_avg(1);
			c = a_A_row_mag_2_avg(2);
		}
		if (a_d == 2)
		{
			a = a_A_row_mag_3_avg(0);
			b = a_A_row_mag_3_avg(1);
			c = a_A_row_mag_3_avg(2);
		}

		W_normalized_sph(0) = W_Sph(0);
		W_normalized_sph(1) = W_Sph(1)*a;
		W_normalized_sph(2) = W_Sph(2)*b;
		W_normalized_sph(3) = W_Sph(3)*c;
		W_normalized_sph(4) = W_Sph(4);
		W_normalized_sph(5) = W_Sph(5)*a;
		W_normalized_sph(6) = W_Sph(6)*b;
		W_normalized_sph(7) = W_Sph(7)*c;
		
	}
	PROTO_KERNEL_END(W_Sph_to_W_normalized_sph_calcF, W_Sph_to_W_normalized_sph_calc)


	void W_Sph_to_W_normalized_sph(BoxData<double,NUMCOMPS>& W_normalized_sph,
	                    const BoxData<double,NUMCOMPS>& W_Sph,
						BoxData<double,DIM>& a_A_row_mag_1_avg,
						BoxData<double,DIM>& a_A_row_mag_2_avg,
						BoxData<double,DIM>& a_A_row_mag_3_avg,
	                    int a_d)
	{
		forallInPlace_p(W_Sph_to_W_normalized_sph_calc, W_normalized_sph, W_Sph, a_A_row_mag_1_avg, a_A_row_mag_2_avg, a_A_row_mag_3_avg, a_d);
	}



	void Regular_map_filling_func(MHDLevelDataState& a_state){

		HDF5Handler h5;
		// for (auto dit : a_state.m_Jacobian_ave){		
		// 	MHD_Mapping::Jacobian_Ave_calc((a_state.m_Jacobian_ave)[dit],a_state.m_dx,a_state.m_dy,a_state.m_dz,a_state.m_U[dit].box().grow(1));
		// 	MHD_Mapping::N_ave_f_calc_func((a_state.m_N_ave_f)[dit],a_state.m_dx,a_state.m_dy,a_state.m_dz);			
		// }
		// (a_state.m_Jacobian_ave).exchange();

		double a_dx = a_state.m_dx;
		double a_dy = a_state.m_dy;
		double a_dz = a_state.m_dz;
		double dxd[3] = {a_dx, a_dy, a_dz};
		for (auto dit : a_state.m_X_corners){
			Box bxmap = a_state.m_X_corners[dit].box();
        	for (int dir = 0;dir <DIM;dir++)
        	{
            	bxmap = bxmap.extrude(dir);
        	}
			BoxData<double,DIM> eta(bxmap);
			MHD_Mapping::etaCorner_calc(eta,bxmap,a_dx, a_dy, a_dz);
			BoxData<double,DIM> X(bxmap);		
			MHD_Mapping::eta_to_x_calc(X, eta, bxmap);

			Array<BoxData<double,DIM>,DIM> NT;
			
			for (int dir = 0; dir < DIM; dir++)
			{
				
				NT[dir] = Operator::cofactor(X,dir);
				NT[dir].copyTo(a_state.m_NT[dir][dit]);
			}
			BoxData<double> J;
			{
				J = Operator::jacobian(X,NT);
				J.copyTo(a_state.m_J[dit]);
				#if DIM==2
					a_state.m_J[dit] *= 1.0/(a_dx*a_dy);
				#endif
				#if DIM==3
					a_state.m_J[dit] *= 1.0/(a_dx*a_dy*a_dz);
				#endif
				// h5.writePatch(1,J,"J");
			}

		}
	}


	PROTO_KERNEL_START
	void Nineto33_calcF( Point& a_pt,
	                          Var<double,DIM,MEM,DIM>& a_A_face_avg,
							  Var<double,DIM*DIM>& a_A_avg)
	{
		a_A_face_avg(0,0) = a_A_avg(0);
		a_A_face_avg(0,1) = a_A_avg(1);
		a_A_face_avg(0,2) = a_A_avg(2);
		a_A_face_avg(1,0) = a_A_avg(3);
		a_A_face_avg(1,1) = a_A_avg(4);
		a_A_face_avg(1,2) = a_A_avg(5);
		a_A_face_avg(2,0) = a_A_avg(6);
		a_A_face_avg(2,1) = a_A_avg(7);
		a_A_face_avg(2,2) = a_A_avg(8);
	}
	PROTO_KERNEL_END(Nineto33_calcF, Nineto33_calc)

	void Nineto33(BoxData<double,DIM,MEM,DIM>& a_A_face_avg,
						BoxData<double,DIM*DIM>& a_A_avg)
	{
		forallInPlace_p(Nineto33_calc, a_A_face_avg, a_A_avg);
	}


	void Spherical_map_filling_func(MHDLevelDataState& a_state)
	{
		bool exchanged_yet = false;
		bool r_dir_turn = false;
		#if DIM == 3
		for (auto dit : a_state.m_detAA_avg){	
			MHD_Mapping::Spherical_map_calc_func((a_state.m_Jacobian_ave)[dit], (a_state.m_A_avg)[dit], (a_state.m_A_inv_avg)[dit], (a_state.m_A_1_avg)[dit], (a_state.m_A_2_avg)[dit], (a_state.m_A_3_avg)[dit], (a_state.m_A_inv_1_avg)[dit], (a_state.m_A_inv_2_avg)[dit], (a_state.m_A_inv_3_avg)[dit], (a_state.m_detAA_avg)[dit], (a_state.m_detAA_inv_avg)[dit], (a_state.m_r2rdot_avg)[dit], (a_state.m_detA_avg)[dit], (a_state.m_A_row_mag_avg)[dit], (a_state.m_r2detA_1_avg)[dit], (a_state.m_r2detAA_1_avg)[dit], (a_state.m_r2detAn_1_avg)[dit], (a_state.m_n_1_avg)[dit], (a_state.m_A_row_mag_1_avg)[dit], (a_state.m_rrdotdetA_2_avg)[dit], (a_state.m_rrdotdetAA_2_avg)[dit], (a_state.m_rrdotd3ncn_2_avg)[dit], (a_state.m_A_row_mag_2_avg)[dit], (a_state.m_rrdotdetA_3_avg)[dit], (a_state.m_rrdotdetAA_3_avg)[dit], (a_state.m_rrdotncd2n_3_avg)[dit], (a_state.m_A_row_mag_3_avg)[dit],a_state.m_dx,a_state.m_dy,a_state.m_dz, exchanged_yet, r_dir_turn);
		}	

		if (inputs.grid_type_global == 2){
			(a_state.m_Jacobian_ave).exchange();
			(a_state.m_A_1_avg).exchange();
			(a_state.m_A_2_avg).exchange();
			(a_state.m_A_3_avg).exchange();
			(a_state.m_A_inv_1_avg).exchange();
			(a_state.m_A_inv_2_avg).exchange();
			(a_state.m_A_inv_3_avg).exchange();
			(a_state.m_detAA_avg).exchange();
			(a_state.m_detAA_inv_avg).exchange();
			(a_state.m_r2rdot_avg).exchange();
			(a_state.m_detA_avg).exchange();
			(a_state.m_r2detA_1_avg).exchange();
			(a_state.m_r2detAA_1_avg).exchange();
			(a_state.m_r2detAn_1_avg).exchange();
			(a_state.m_rrdotdetA_2_avg).exchange();
			(a_state.m_rrdotdetAA_2_avg).exchange();
			(a_state.m_rrdotd3ncn_2_avg).exchange();
			(a_state.m_rrdotdetA_3_avg).exchange();
			(a_state.m_rrdotdetAA_3_avg).exchange();
			(a_state.m_rrdotncd2n_3_avg).exchange();
			exchanged_yet = true;
		}
		for (auto dit : a_state.m_detAA_avg){		
			MHD_Mapping::Spherical_map_calc_func((a_state.m_Jacobian_ave)[dit], (a_state.m_A_avg)[dit], (a_state.m_A_inv_avg)[dit], (a_state.m_A_1_avg)[dit], (a_state.m_A_2_avg)[dit], (a_state.m_A_3_avg)[dit], (a_state.m_A_inv_1_avg)[dit], (a_state.m_A_inv_2_avg)[dit], (a_state.m_A_inv_3_avg)[dit], (a_state.m_detAA_avg)[dit], (a_state.m_detAA_inv_avg)[dit], (a_state.m_r2rdot_avg)[dit], (a_state.m_detA_avg)[dit], (a_state.m_A_row_mag_avg)[dit], (a_state.m_r2detA_1_avg)[dit], (a_state.m_r2detAA_1_avg)[dit], (a_state.m_r2detAn_1_avg)[dit], (a_state.m_n_1_avg)[dit], (a_state.m_A_row_mag_1_avg)[dit], (a_state.m_rrdotdetA_2_avg)[dit], (a_state.m_rrdotdetAA_2_avg)[dit], (a_state.m_rrdotd3ncn_2_avg)[dit], (a_state.m_A_row_mag_2_avg)[dit], (a_state.m_rrdotdetA_3_avg)[dit], (a_state.m_rrdotdetAA_3_avg)[dit], (a_state.m_rrdotncd2n_3_avg)[dit], (a_state.m_A_row_mag_3_avg)[dit],a_state.m_dx,a_state.m_dy,a_state.m_dz, exchanged_yet, r_dir_turn);
		}
		exchanged_yet = false;
		r_dir_turn = true;

		for (auto dit : a_state.m_detAA_avg){		
			MHD_Mapping::Spherical_map_calc_func((a_state.m_Jacobian_ave)[dit], (a_state.m_A_avg)[dit], (a_state.m_A_inv_avg)[dit], (a_state.m_A_1_avg)[dit], (a_state.m_A_2_avg)[dit], (a_state.m_A_3_avg)[dit], (a_state.m_A_inv_1_avg)[dit], (a_state.m_A_inv_2_avg)[dit], (a_state.m_A_inv_3_avg)[dit], (a_state.m_detAA_avg)[dit], (a_state.m_detAA_inv_avg)[dit], (a_state.m_r2rdot_avg)[dit], (a_state.m_detA_avg)[dit], (a_state.m_A_row_mag_avg)[dit], (a_state.m_r2detA_1_avg)[dit], (a_state.m_r2detAA_1_avg)[dit], (a_state.m_r2detAn_1_avg)[dit], (a_state.m_n_1_avg)[dit], (a_state.m_A_row_mag_1_avg)[dit], (a_state.m_rrdotdetA_2_avg)[dit], (a_state.m_rrdotdetAA_2_avg)[dit], (a_state.m_rrdotd3ncn_2_avg)[dit], (a_state.m_A_row_mag_2_avg)[dit], (a_state.m_rrdotdetA_3_avg)[dit], (a_state.m_rrdotdetAA_3_avg)[dit], (a_state.m_rrdotncd2n_3_avg)[dit], (a_state.m_A_row_mag_3_avg)[dit],a_state.m_dx,a_state.m_dy,a_state.m_dz, exchanged_yet, r_dir_turn);
		}

		

		#endif
	}


	void Spherical_map_filling_func2(MHDLevelDataState& a_state)
	{
		//Filling data for Phil's operators
		for (auto dit : a_state.m_A_1_avg){		
			a_state.m_A_1_avg[dit].copyTo(a_state.m_A_face_avg[0][dit]);
			a_state.m_A_2_avg[dit].copyTo(a_state.m_A_face_avg[1][dit]);
			a_state.m_A_3_avg[dit].copyTo(a_state.m_A_face_avg[2][dit]);
			a_state.m_r2detA_1_avg[dit].copyTo(a_state.m_Dr_detA_avg[0][dit]);
			a_state.m_rrdotdetA_2_avg[dit].copyTo(a_state.m_Dr_detA_avg[1][dit]);
			a_state.m_rrdotdetA_3_avg[dit].copyTo(a_state.m_Dr_detA_avg[2][dit]);
			a_state.m_r2detAn_1_avg[dit].copyTo(a_state.m_Dr_adjA_avg[0][dit]);
			a_state.m_rrdotd3ncn_2_avg[dit].copyTo(a_state.m_Dr_adjA_avg[1][dit]);
			a_state.m_rrdotncd2n_3_avg[dit].copyTo(a_state.m_Dr_adjA_avg[2][dit]);

		}
	}

	











	PROTO_KERNEL_START
	void Cartesian_to_Spherical_calcF( const Point& a_pt,
	                          Var<double,NUMCOMPS>& W_sph,
	                          Var<double,NUMCOMPS>& W_cart,
							  const Var<double,DIM>& a_x_sph)
	{
		
		double r = a_x_sph(0);
		double theta = a_x_sph(1);
		double phi = a_x_sph(2);

		double A0 = sin(theta)*cos(phi);
		double A1 = sin(theta)*sin(phi);
		double A2 = cos(theta);
		double A3 = cos(theta)*cos(phi);
		double A4 = cos(theta)*sin(phi);
		double A5 = -sin(theta);
		double A6 = -sin(phi);
		double A7 = cos(phi);
		double A8 = 0.0;

		
		W_sph(0) = W_cart(0);
		W_sph(4) = W_cart(4);
		W_sph(1) = A0*W_cart(1) + A1*W_cart(2) + A2*W_cart(3);
		W_sph(2) = A3*W_cart(1) + A4*W_cart(2) + A5*W_cart(3);
		W_sph(3) = A6*W_cart(1) + A7*W_cart(2) + A8*W_cart(3);
		W_sph(5) = A0*W_cart(5) + A1*W_cart(6) + A2*W_cart(7);
		W_sph(6) = A3*W_cart(5) + A4*W_cart(6) + A5*W_cart(7);
		W_sph(7) = A6*W_cart(5) + A7*W_cart(6) + A8*W_cart(7);
	}
	PROTO_KERNEL_END(Cartesian_to_Spherical_calcF, Cartesian_to_Spherical_calc)

	void Cartesian_to_Spherical(BoxData<double,NUMCOMPS>& a_W_Sph,
	                    		const BoxData<double,NUMCOMPS>& a_W_Cart,
	                    		const BoxData<double,DIM>& a_x_Sph)
	{
		forallInPlace_p(Cartesian_to_Spherical_calc, a_W_Sph, a_W_Cart, a_x_Sph);
	}	


	PROTO_KERNEL_START
	void Spherical_to_Cartesian_calcF( const Point& a_pt,
	                          Var<double,NUMCOMPS>& W_cart,
	                          Var<double,NUMCOMPS>& W_sph,
							  const Var<double,DIM>& a_x_sph)
	{
		
		double r = a_x_sph(0);
		double theta = a_x_sph(1);
		double phi = a_x_sph(2);

		double A0 = sin(theta)*cos(phi);
		double A1 = cos(theta)*cos(phi);
		double A2 = -sin(phi);
		double A3 = sin(theta)*sin(phi);
		double A4 = cos(theta)*sin(phi);
		double A5 = cos(phi);
		double A6 = cos(theta);
		double A7 = -sin(theta);
		double A8 = 0.0;

		
		W_cart(0) = W_sph(0);
		W_cart(4) = W_sph(4);
		W_cart(1) = A0*W_sph(1) + A1*W_sph(2) + A2*W_sph(3);
		W_cart(2) = A3*W_sph(1) + A4*W_sph(2) + A5*W_sph(3);
		W_cart(3) = A6*W_sph(1) + A7*W_sph(2) + A8*W_sph(3);
		W_cart(5) = A0*W_sph(5) + A1*W_sph(6) + A2*W_sph(7);
		W_cart(6) = A3*W_sph(5) + A4*W_sph(6) + A5*W_sph(7);
		W_cart(7) = A6*W_sph(5) + A7*W_sph(6) + A8*W_sph(7);
	}
	PROTO_KERNEL_END(Spherical_to_Cartesian_calcF, Spherical_to_Cartesian_calc)

	void Spherical_to_Cartesian(BoxData<double,NUMCOMPS>& a_W_Cart,
	                    		const BoxData<double,NUMCOMPS>& a_W_Sph,
	                    		const BoxData<double,DIM>& a_x_Sph)
	{
		forallInPlace_p(Spherical_to_Cartesian_calc, a_W_Cart, a_W_Sph, a_x_Sph);
	}



	PROTO_KERNEL_START
	void get_sph_coords_fc_calcF(const Point& a_pt,
							Var<double,DIM>& a_x_sph,
	                        const double a_dx,
	                   		const double a_dy,
	                   		const double a_dz,
	                   		int a_d)
	{
		double eta[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0) dxd = a_dx;
			if (i == 1) dxd = a_dy;
			if (i == 2) dxd = a_dz;
			if (i == a_d) {
				eta[i] = a_pt[i]*dxd;
			} else {
				eta[i] = a_pt[i]*dxd + 0.5*dxd;
			}
		}

		double R_t = (inputs.r_out*c_AU - inputs.r_in*c_AU)/(exp(inputs.C_rad) - 1.0);

		double r = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*eta[0]) - 1.0);
		double theta = c_PI*eta[1];
		double phi = 2.0*c_PI*eta[2];

		if (theta < 0){
			theta = -theta;
			if (phi < c_PI){
				phi = phi + c_PI;
			} else {
				phi = phi - c_PI;
			}
		}
		if (theta > c_PI){
			theta = c_PI - (theta - c_PI);
			if (phi < c_PI){
				phi = phi + c_PI;
			} else {
				phi = phi - c_PI;
			}
		}
		a_x_sph(0) = r;  //r
		a_x_sph(1) = theta; //theta
		a_x_sph(2) = phi; //phi
	}
	PROTO_KERNEL_END(get_sph_coords_fc_calcF, get_sph_coords_fc_calc)


	void get_sph_coords_fc(BoxData<double,DIM>& a_x_sph,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz,
	                    int a_d)
	{
		forallInPlace_p(get_sph_coords_fc_calc, a_bx, a_x_sph, a_dx, a_dy, a_dz, a_d);
	}

	

	PROTO_KERNEL_START
	void get_sph_coords_cc_calcF(const Point& a_pt,
							Var<double,DIM>& a_x_sph,
	                        const double a_dx,
	                   		const double a_dy,
	                   		const double a_dz)
	{
		double R_t = (inputs.r_out*c_AU - inputs.r_in*c_AU)/(exp(inputs.C_rad) - 1.0);
		double eta_here[3];
		double eta_ahead[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0) dxd = a_dx;
			if (i == 1) dxd = a_dy;
			if (i == 2) dxd = a_dz;
			eta_here[i] = a_pt[i]*dxd;
			eta_ahead[i] = a_pt[i]*dxd + dxd;
		}

		double r_here = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*eta_here[0]) - 1.0);
		double theta_here = c_PI*eta_here[1];
		double phi_here = 2.0*c_PI*eta_here[2];

		double r_ahead = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*eta_ahead[0]) - 1.0);
		double theta_ahead = c_PI*eta_ahead[1];
		double phi_ahead = 2.0*c_PI*eta_ahead[2];

		double r = 0.5*(r_here + r_ahead);
		double theta = 0.5*(theta_here + theta_ahead);
		double phi = 0.5*(phi_here + phi_ahead);

		if (theta < 0){
			theta = -theta;
			if (phi < c_PI){
				phi = phi + c_PI;
			} else {
				phi = phi - c_PI;
			}
		}
		if (theta > c_PI){
			theta = c_PI - (theta - c_PI);
			if (phi < c_PI){
				phi = phi + c_PI;
			} else {
				phi = phi - c_PI;
			}
		}
		a_x_sph(0) = r;  //r
		a_x_sph(1) = theta; //theta
		a_x_sph(2) = phi; //phi
	}
	PROTO_KERNEL_END(get_sph_coords_cc_calcF, get_sph_coords_cc_calc)


	void get_sph_coords_cc(BoxData<double,DIM>& a_x_sph,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz)
	{
		forallInPlace_p(get_sph_coords_cc_calc, a_bx, a_x_sph, a_dx, a_dy, a_dz);
	}
	


	PROTO_KERNEL_START
	void get_cell_volume_calcF(const Point& a_pt,
							Var<double,1>& a_V,
	                        const double a_dx,
	                   		const double a_dy,
	                   		const double a_dz)
	{
		double R_t = (inputs.r_out*c_AU - inputs.r_in*c_AU)/(exp(inputs.C_rad) - 1.0);
		double eta_here[3];
		double eta_ahead[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0) dxd = a_dx;
			if (i == 1) dxd = a_dy;
			if (i == 2) dxd = a_dz;
			eta_here[i] = a_pt[i]*dxd;
			eta_ahead[i] = a_pt[i]*dxd + dxd;
		}

		double r_here = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*eta_here[0]) - 1.0);
		double theta_here = c_PI*eta_here[1];
		double phi_here = 2.0*c_PI*eta_here[2];

		double r_ahead = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*eta_ahead[0]) - 1.0);
		double theta_ahead = c_PI*eta_ahead[1];
		double phi_ahead = 2.0*c_PI*eta_ahead[2];

		double volume = (1.0/3.0)*(pow(r_ahead,3)-pow(r_here,3))*(cos(theta_here)-cos(theta_ahead))*(phi_ahead-phi_here);
		
		a_V(0) = abs(volume);  
		
	}
	PROTO_KERNEL_END(get_cell_volume_calcF, get_cell_volume_calc)


	void get_cell_volume(BoxData<double,1>& a_V,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz)
	{
		forallInPlace_p(get_cell_volume_calc, a_bx, a_V, a_dx, a_dy, a_dz);
	}


	PROTO_KERNEL_START
	void get_face_area_calcF(const Point& a_pt,
							Var<double,DIM>& a_A,
	                        const double a_dx,
	                   		const double a_dy,
	                   		const double a_dz)
	{
		double R_t = (inputs.r_out*c_AU - inputs.r_in*c_AU)/(exp(inputs.C_rad) - 1.0);
		double eta_here[3];
		double eta_ahead[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0) dxd = a_dx;
			if (i == 1) dxd = a_dy;
			if (i == 2) dxd = a_dz;
			eta_here[i] = a_pt[i]*dxd;
			eta_ahead[i] = a_pt[i]*dxd + dxd;
		}

		double r_here = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*eta_here[0]) - 1.0);
		double theta_here = c_PI*eta_here[1];
		double phi_here = 2.0*c_PI*eta_here[2];

		double r_ahead = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*eta_ahead[0]) - 1.0);
		double theta_ahead = c_PI*eta_ahead[1];
		double phi_ahead = 2.0*c_PI*eta_ahead[2];

		double A_r = pow(r_here,2)*(cos(theta_here)-cos(theta_ahead))*(phi_ahead-phi_here);
		double A_theta = 0.5*sin(theta_here)*(pow(r_ahead,2) - pow(r_here,2))*(phi_ahead-phi_here);
		double A_phi = 0.5*(pow(r_ahead,2) - pow(r_here,2))*(theta_ahead-theta_here);
		
		a_A(0) = abs(A_r);  
		a_A(1) = abs(A_theta);  
		a_A(2) = abs(A_phi);  
		
	}
	PROTO_KERNEL_END(get_face_area_calcF, get_face_area_calc)


	void get_face_area(BoxData<double,DIM>& a_A,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz)
	{
		forallInPlace_p(get_face_area_calc, a_bx, a_A, a_dx, a_dy, a_dz);
	}

	PROTO_KERNEL_START
	void get_delta_sph_coords_calcF(const Point& a_pt,
							Var<double,DIM>& a_dx_sph,
	                        const double a_dx,
	                   		const double a_dy,
	                   		const double a_dz)
	{
		double R_t = (inputs.r_out*c_AU - inputs.r_in*c_AU)/(exp(inputs.C_rad) - 1.0);
		double eta_here[3];
		double eta_ahead[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0) dxd = a_dx;
			if (i == 1) dxd = a_dy;
			if (i == 2) dxd = a_dz;
			eta_here[i] = a_pt[i]*dxd;
			eta_ahead[i] = a_pt[i]*dxd + dxd;
		}

		double r_here = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*eta_here[0]) - 1.0);
		double theta_here = c_PI*eta_here[1];
		double phi_here = 2.0*c_PI*eta_here[2];

		double r_ahead = inputs.r_in*c_AU + R_t*(exp(inputs.C_rad*eta_ahead[0]) - 1.0);
		double theta_ahead = c_PI*eta_ahead[1];
		double phi_ahead = 2.0*c_PI*eta_ahead[2];

		double dr = r_ahead - r_here;
		double dtheta = theta_ahead - theta_here;
		double dphi = phi_ahead - phi_here;
		
		a_dx_sph(0) = dr;  //dr
		a_dx_sph(1) = dtheta; //dtheta
		a_dx_sph(2) = dphi; //dphi
	}
	PROTO_KERNEL_END(get_delta_sph_coords_calcF, get_delta_sph_coords_calc)


	void get_delta_sph_coords(BoxData<double,DIM>& a_dx_sph,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz)
	{
		forallInPlace_p(get_delta_sph_coords_calc, a_bx, a_dx_sph, a_dx, a_dy, a_dz);
	}
	PROTO_KERNEL_START
	void Correct_V_theta_phi_at_poles_calcF( const Point& a_pt,
	                          Var<double,NUMCOMPS>& a_U_Sph_ave,
	                          const double a_dx,
	                          const double a_dy,
	                          const double a_dz)
	{

		double E1, E2, E3;
		E1 = (a_pt[0] + 0.5)*a_dx;
		E2 = (a_pt[1] + 0.5)*a_dy;
		E3 = (a_pt[2] + 0.5)*a_dz;

		if (E2 < 0.0 || E2 > 1.0){
			a_U_Sph_ave(2) *= -1.0;
			a_U_Sph_ave(3) *= -1.0;

			a_U_Sph_ave(6) *= -1.0;
			a_U_Sph_ave(7) *= -1.0;
		}

	}
	PROTO_KERNEL_END(Correct_V_theta_phi_at_poles_calcF, Correct_V_theta_phi_at_poles_calc)

	void Correct_V_theta_phi_at_poles(BoxData<double,NUMCOMPS>& a_U_Sph_ave,
	                             const double a_dx,
	                             const double a_dy,
	                             const double a_dz)
	{
		forallInPlace_p(Correct_V_theta_phi_at_poles_calc, a_U_Sph_ave, a_dx, a_dy, a_dz);
	}

	void Spherical_2O_map_filling_func(MHDLevelDataState& a_state)
	{
		for (auto dit : a_state.m_cell_volume){	
			Box dbx0 = a_state.m_cell_volume[dit].box();
			MHD_Mapping::get_cell_volume(a_state.m_cell_volume[dit],dbx0,a_state.m_dx,a_state.m_dy,a_state.m_dz);
			MHD_Mapping::get_face_area(a_state.m_face_area[dit],dbx0,a_state.m_dx,a_state.m_dy,a_state.m_dz);
			MHD_Mapping::get_delta_sph_coords(a_state.m_dx_sph[dit],dbx0,a_state.m_dx,a_state.m_dy,a_state.m_dz);
			MHD_Mapping::get_sph_coords_cc(a_state.m_x_sph_cc[dit],dbx0,a_state.m_dx, a_state.m_dy, a_state.m_dz);
			MHD_Mapping::get_sph_coords_fc(a_state.m_x_sph_fc_1[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz, 0);
			MHD_Mapping::get_sph_coords_fc(a_state.m_x_sph_fc_2[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz, 1);
			MHD_Mapping::get_sph_coords_fc(a_state.m_x_sph_fc_3[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz, 2);
		}
	}
}
