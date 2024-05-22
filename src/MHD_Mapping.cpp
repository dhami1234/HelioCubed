#include "Proto.H"
#include "MHD_Mapping.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
// #include "Proto_WriteBoxData.H"
#include "MHDOp.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Constants.H"
#include "MHD_Operator.H"
// #include "MHDLevelDataRK4.H"
extern Parsefrominputs inputs;

// constexpr MemType MEM = MEMTYPE_DEFAULT;

typedef BoxData<double, 1, HOST> Scalar;
typedef BoxData<double, NUMCOMPS, HOST> Vector;
/// @brief MHD_Mapping namespace
namespace MHD_Mapping
{

	PROTO_KERNEL_START
	void out_data_joinF(Var<double, NUMCOMPS + DIM> &a_out_data,
						const Var<double, DIM> &a_phys_coords,
						const Var<double, NUMCOMPS> &a_W)
	{
		for (int i = 0; i < DIM; i++)
		{
			a_out_data(i) = a_phys_coords(i) / c_AU; // AU
		}
		a_out_data(DIM + 0) = a_W(0) / c_MP; // /cm^3
		a_out_data(DIM + 1) = a_W(1) / 1e5;	 // km/s
		a_out_data(DIM + 2) = a_W(2) / 1e5;	 // km/s
		a_out_data(DIM + 3) = a_W(3) / 1e5;	 // km/s
		a_out_data(DIM + 4) = a_W(4);		 // dyne/cm*2
		a_out_data(DIM + 5) = a_W(5);		 // Gauss
		a_out_data(DIM + 6) = a_W(6);		 // Gauss
		a_out_data(DIM + 7) = a_W(7);		 // Gauss
		#if TURB == 1
		a_out_data(DIM + iZ2) = a_W(iZ2);
		a_out_data(DIM + iSIGMA) = a_W(iSIGMA);
		a_out_data(DIM + iLAMBDA) = a_W(iLAMBDA);
		#endif
		
	}
	PROTO_KERNEL_END(out_data_joinF, out_data_join)

	void out_data_calc(BoxData<double, NUMCOMPS + DIM> &a_out_data,
					   const BoxData<double, DIM> &a_phys_coords,
					   const BoxData<double, NUMCOMPS> &a_W)
	{
		a_out_data = forall<double, NUMCOMPS + DIM>(out_data_join, a_phys_coords, a_W);
	}


	PROTO_KERNEL_START
	void Cartesian_to_Spherical_calcF(const Point &a_pt,
									  Var<double, NUMCOMPS> &W_sph,
									  Var<double, NUMCOMPS> &W_cart,
									  const Var<double, DIM> &a_x_sph)
	{

		double r = a_x_sph(0);
		double theta = a_x_sph(1);
		double phi = a_x_sph(2);

		double A0 = sin(theta) * cos(phi);
		double A1 = sin(theta) * sin(phi);
		double A2 = cos(theta);
		double A3 = cos(theta) * cos(phi);
		double A4 = cos(theta) * sin(phi);
		double A5 = -sin(theta);
		double A6 = -sin(phi);
		double A7 = cos(phi);
		double A8 = 0.0;

		W_sph(0) = W_cart(0);
		W_sph(4) = W_cart(4);
		W_sph(1) = A0 * W_cart(1) + A1 * W_cart(2) + A2 * W_cart(3);
		W_sph(2) = A3 * W_cart(1) + A4 * W_cart(2) + A5 * W_cart(3);
		W_sph(3) = A6 * W_cart(1) + A7 * W_cart(2) + A8 * W_cart(3);
		W_sph(5) = A0 * W_cart(5) + A1 * W_cart(6) + A2 * W_cart(7);
		W_sph(6) = A3 * W_cart(5) + A4 * W_cart(6) + A5 * W_cart(7);
		W_sph(7) = A6 * W_cart(5) + A7 * W_cart(6) + A8 * W_cart(7);

		#if TURB == 1
		W_sph(iZ2) = W_cart(iZ2);
		W_sph(iSIGMA) = W_cart(iSIGMA);
		W_sph(iLAMBDA) = W_cart(iLAMBDA);
		#endif
	}
	PROTO_KERNEL_END(Cartesian_to_Spherical_calcF, Cartesian_to_Spherical_calc)

	void Cartesian_to_Spherical(BoxData<double, NUMCOMPS> &a_W_Sph,
								const BoxData<double, NUMCOMPS> &a_W_Cart,
								const BoxData<double, DIM> &a_x_Sph)
	{
		forallInPlace_p(Cartesian_to_Spherical_calc, a_W_Sph, a_W_Cart, a_x_Sph);
	}

	PROTO_KERNEL_START
	void Spherical_to_Cartesian_calcF(const Point &a_pt,
									  Var<double, NUMCOMPS> &W_cart,
									  Var<double, NUMCOMPS> &W_sph,
									  const Var<double, DIM> &a_x_sph)
	{

		double r = a_x_sph(0);
		double theta = a_x_sph(1);
		double phi = a_x_sph(2);

		double A0 = sin(theta) * cos(phi);
		double A1 = cos(theta) * cos(phi);
		double A2 = -sin(phi);
		double A3 = sin(theta) * sin(phi);
		double A4 = cos(theta) * sin(phi);
		double A5 = cos(phi);
		double A6 = cos(theta);
		double A7 = -sin(theta);
		double A8 = 0.0;

		W_cart(0) = W_sph(0);
		W_cart(4) = W_sph(4);
		W_cart(1) = A0 * W_sph(1) + A1 * W_sph(2) + A2 * W_sph(3);
		W_cart(2) = A3 * W_sph(1) + A4 * W_sph(2) + A5 * W_sph(3);
		W_cart(3) = A6 * W_sph(1) + A7 * W_sph(2) + A8 * W_sph(3);
		W_cart(5) = A0 * W_sph(5) + A1 * W_sph(6) + A2 * W_sph(7);
		W_cart(6) = A3 * W_sph(5) + A4 * W_sph(6) + A5 * W_sph(7);
		W_cart(7) = A6 * W_sph(5) + A7 * W_sph(6) + A8 * W_sph(7);

		#if TURB == 1
		W_cart(iZ2) = W_sph(iZ2);
		W_cart(iSIGMA) = W_sph(iSIGMA);
		W_cart(iLAMBDA) = W_sph(iLAMBDA);
		#endif
	}
	PROTO_KERNEL_END(Spherical_to_Cartesian_calcF, Spherical_to_Cartesian_calc)

	void Spherical_to_Cartesian(BoxData<double, NUMCOMPS> &a_W_Cart,
								const BoxData<double, NUMCOMPS> &a_W_Sph,
								const BoxData<double, DIM> &a_x_Sph)
	{
		forallInPlace_p(Spherical_to_Cartesian_calc, a_W_Cart, a_W_Sph, a_x_Sph);
	}

	PROTO_KERNEL_START
	void get_sph_coords_fc_calcF(const Point &a_pt,
								 Var<double, DIM> &a_x_sph,
								 const double a_dx,
								 const double a_dy,
								 const double a_dz,
								 int a_d)
	{
		double eta[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0)
				dxd = a_dx;
			if (i == 1)
				dxd = a_dy;
			if (i == 2)
				dxd = a_dz;
			if (i == a_d)
			{
				eta[i] = a_pt[i] * dxd;
			}
			else
			{
				eta[i] = a_pt[i] * dxd + 0.5 * dxd;
			}
		}

		double R_t = (inputs.r_out * c_AU - inputs.r_in * c_AU) / (exp(inputs.C_rad) - 1.0);

		double r = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta[0]) - 1.0);
		double theta = c_PI * eta[1];
		double phi = 2.0 * c_PI * eta[2];

		if (theta < 0)
		{
			theta = -theta;
			if (phi < c_PI)
			{
				phi = phi + c_PI;
			}
			else
			{
				phi = phi - c_PI;
			}
		}
		if (theta > c_PI)
		{
			theta = c_PI - (theta - c_PI);
			if (phi < c_PI)
			{
				phi = phi + c_PI;
			}
			else
			{
				phi = phi - c_PI;
			}
		}
		a_x_sph(0) = r;		// r
		a_x_sph(1) = theta; // theta
		a_x_sph(2) = phi;	// phi
	}
	PROTO_KERNEL_END(get_sph_coords_fc_calcF, get_sph_coords_fc_calc)

	void get_sph_coords_fc(BoxData<double, DIM> &a_x_sph,
						   const Box &a_bx,
						   const double a_dx,
						   const double a_dy,
						   const double a_dz,
						   int a_d)
	{
		forallInPlace_p(get_sph_coords_fc_calc, a_bx, a_x_sph, a_dx, a_dy, a_dz, a_d);
	}

	PROTO_KERNEL_START
	void get_cart_coords_fc_calcF(const Point &a_pt,
								 Var<double, DIM> &a_x_cart,
								 const double a_dx,
								 const double a_dy,
								 const double a_dz,
								 int a_d)
	{
		double eta[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0)
				dxd = a_dx;
			if (i == 1)
				dxd = a_dy;
			if (i == 2)
				dxd = a_dz;
			if (i == a_d)
			{
				eta[i] = a_pt[i] * dxd;
			}
			else
			{
				eta[i] = a_pt[i] * dxd + 0.5 * dxd;
			}
		}

		double R_t = (inputs.r_out * c_AU - inputs.r_in * c_AU) / (exp(inputs.C_rad) - 1.0);

		double r = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta[0]) - 1.0);
		double theta = c_PI * eta[1];
		double phi = 2.0 * c_PI * eta[2];

		if (theta < 0)
		{
			theta = -theta;
			if (phi < c_PI)
			{
				phi = phi + c_PI;
			}
			else
			{
				phi = phi - c_PI;
			}
		}
		if (theta > c_PI)
		{
			theta = c_PI - (theta - c_PI);
			if (phi < c_PI)
			{
				phi = phi + c_PI;
			}
			else
			{
				phi = phi - c_PI;
			}
		}
		a_x_cart(0) = r*sin(theta)*cos(phi);	// x
		a_x_cart(1) = r*sin(theta)*sin(phi);    // y
		a_x_cart(2) = r*cos(theta);				// z
	}
	PROTO_KERNEL_END(get_cart_coords_fc_calcF, get_cart_coords_fc_calc)

	void get_cart_coords_fc(BoxData<double, DIM> &a_x_cart,
						   const Box &a_bx,
						   const double a_dx,
						   const double a_dy,
						   const double a_dz,
						   int a_d)
	{
		forallInPlace_p(get_cart_coords_fc_calc, a_bx, a_x_cart, a_dx, a_dy, a_dz, a_d);
	}

	PROTO_KERNEL_START
	void get_sph_coords_cc_calcF(const Point &a_pt,
								 Var<double, DIM> &a_x_sph,
								 const double a_dx,
								 const double a_dy,
								 const double a_dz)
	{
		double R_t = (inputs.r_out * c_AU - inputs.r_in * c_AU) / (exp(inputs.C_rad) - 1.0);
		double eta_here[3];
		double eta_ahead[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0)
				dxd = a_dx;
			if (i == 1)
				dxd = a_dy;
			if (i == 2)
				dxd = a_dz;
			eta_here[i] = a_pt[i] * dxd;
			eta_ahead[i] = a_pt[i] * dxd + dxd;
		}

		double r_here = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_here[0]) - 1.0);
		double theta_here = c_PI * eta_here[1];
		double phi_here = 2.0 * c_PI * eta_here[2];

		double r_ahead = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_ahead[0]) - 1.0);
		double theta_ahead = c_PI * eta_ahead[1];
		double phi_ahead = 2.0 * c_PI * eta_ahead[2];

		double r = 0.5 * (r_here + r_ahead);
		double theta = 0.5 * (theta_here + theta_ahead);
		double phi = 0.5 * (phi_here + phi_ahead);

		if (theta < 0)
		{
			theta = -theta;
			if (phi < c_PI)
			{
				phi = phi + c_PI;
			}
			else
			{
				phi = phi - c_PI;
			}
		}
		if (theta > c_PI)
		{
			theta = c_PI - (theta - c_PI);
			if (phi < c_PI)
			{
				phi = phi + c_PI;
			}
			else
			{
				phi = phi - c_PI;
			}
		}
		a_x_sph(0) = r;		// r
		a_x_sph(1) = theta; // theta
		a_x_sph(2) = phi;	// phi
	}
	PROTO_KERNEL_END(get_sph_coords_cc_calcF, get_sph_coords_cc_calc)

	void get_sph_coords_cc(BoxData<double, DIM> &a_x_sph,
						   const Box &a_bx,
						   const double a_dx,
						   const double a_dy,
						   const double a_dz)
	{
		forallInPlace_p(get_sph_coords_cc_calc, a_bx, a_x_sph, a_dx, a_dy, a_dz);
	}

	PROTO_KERNEL_START
	void get_cart_coords_cc_calcF(const Point &a_pt,
								 Var<double, DIM> &a_x_cart,
								 const double a_dx,
								 const double a_dy,
								 const double a_dz)
	{
		double R_t = (inputs.r_out * c_AU - inputs.r_in * c_AU) / (exp(inputs.C_rad) - 1.0);
		double eta_here[3];
		double eta_ahead[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0)
				dxd = a_dx;
			if (i == 1)
				dxd = a_dy;
			if (i == 2)
				dxd = a_dz;
			eta_here[i] = a_pt[i] * dxd;
			eta_ahead[i] = a_pt[i] * dxd + dxd;
		}

		double r_here = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_here[0]) - 1.0);
		double theta_here = c_PI * eta_here[1];
		double phi_here = 2.0 * c_PI * eta_here[2];

		double r_ahead = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_ahead[0]) - 1.0);
		double theta_ahead = c_PI * eta_ahead[1];
		double phi_ahead = 2.0 * c_PI * eta_ahead[2];

		double r = 0.5 * (r_here + r_ahead);
		double theta = 0.5 * (theta_here + theta_ahead);
		double phi = 0.5 * (phi_here + phi_ahead);

		if (theta < 0)
		{
			theta = -theta;
			if (phi < c_PI)
			{
				phi = phi + c_PI;
			}
			else
			{
				phi = phi - c_PI;
			}
		}
		if (theta > c_PI)
		{
			theta = c_PI - (theta - c_PI);
			if (phi < c_PI)
			{
				phi = phi + c_PI;
			}
			else
			{
				phi = phi - c_PI;
			}
		}
		a_x_cart(0) = r*sin(theta)*cos(phi);	// x
		a_x_cart(1) = r*sin(theta)*sin(phi);    // y
		a_x_cart(2) = r*cos(theta);				// z
	}
	PROTO_KERNEL_END(get_cart_coords_cc_calcF, get_cart_coords_cc_calc)

	void get_cart_coords_cc(BoxData<double, DIM> &a_x_cart,
						   const Box &a_bx,
						   const double a_dx,
						   const double a_dy,
						   const double a_dz)
	{
		forallInPlace_p(get_cart_coords_cc_calc, a_bx, a_x_cart, a_dx, a_dy, a_dz);
	}

	PROTO_KERNEL_START
	void get_cell_volume_calcF(const Point &a_pt,
							   Var<double, 1> &a_V,
							   const double a_dx,
							   const double a_dy,
							   const double a_dz)
	{
		double R_t = (inputs.r_out * c_AU - inputs.r_in * c_AU) / (exp(inputs.C_rad) - 1.0);
		double eta_here[3];
		double eta_ahead[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0)
				dxd = a_dx;
			if (i == 1)
				dxd = a_dy;
			if (i == 2)
				dxd = a_dz;
			eta_here[i] = a_pt[i] * dxd;
			eta_ahead[i] = a_pt[i] * dxd + dxd;
		}

		double r_here = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_here[0]) - 1.0);
		double theta_here = c_PI * eta_here[1];
		double phi_here = 2.0 * c_PI * eta_here[2];

		double r_ahead = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_ahead[0]) - 1.0);
		double theta_ahead = c_PI * eta_ahead[1];
		double phi_ahead = 2.0 * c_PI * eta_ahead[2];

		double volume = (1.0 / 3.0) * (pow(r_ahead, 3) - pow(r_here, 3)) * (cos(theta_here) - cos(theta_ahead)) * (phi_ahead - phi_here);

		a_V(0) = abs(volume);
	}
	PROTO_KERNEL_END(get_cell_volume_calcF, get_cell_volume_calc)

	void get_cell_volume(BoxData<double, 1> &a_V,
						 const Box &a_bx,
						 const double a_dx,
						 const double a_dy,
						 const double a_dz)
	{
		forallInPlace_p(get_cell_volume_calc, a_bx, a_V, a_dx, a_dy, a_dz);
	}

	PROTO_KERNEL_START
	void get_face_area_calcF(const Point &a_pt,
							 Var<double, DIM> &a_A,
							 const double a_dx,
							 const double a_dy,
							 const double a_dz)
	{
		double R_t = (inputs.r_out * c_AU - inputs.r_in * c_AU) / (exp(inputs.C_rad) - 1.0);
		double eta_here[3];
		double eta_ahead[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0)
				dxd = a_dx;
			if (i == 1)
				dxd = a_dy;
			if (i == 2)
				dxd = a_dz;
			eta_here[i] = a_pt[i] * dxd;
			eta_ahead[i] = a_pt[i] * dxd + dxd;
		}

		double r_here = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_here[0]) - 1.0);
		double theta_here = c_PI * eta_here[1];
		double phi_here = 2.0 * c_PI * eta_here[2];

		double r_ahead = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_ahead[0]) - 1.0);
		double theta_ahead = c_PI * eta_ahead[1];
		double phi_ahead = 2.0 * c_PI * eta_ahead[2];

		double A_r = pow(r_here, 2) * (cos(theta_here) - cos(theta_ahead)) * (phi_ahead - phi_here);
		double A_theta = 0.5 * sin(theta_here) * (pow(r_ahead, 2) - pow(r_here, 2)) * (phi_ahead - phi_here);
		double A_phi = 0.5 * (pow(r_ahead, 2) - pow(r_here, 2)) * (theta_ahead - theta_here);

		a_A(0) = abs(A_r);
		a_A(1) = abs(A_theta);
		a_A(2) = abs(A_phi);
	}
	PROTO_KERNEL_END(get_face_area_calcF, get_face_area_calc)

	void get_face_area(BoxData<double, DIM> &a_A,
					   const Box &a_bx,
					   const double a_dx,
					   const double a_dy,
					   const double a_dz)
	{
		forallInPlace_p(get_face_area_calc, a_bx, a_A, a_dx, a_dy, a_dz);
	}

	PROTO_KERNEL_START
	void get_delta_sph_coords_calcF(const Point &a_pt,
									Var<double, DIM> &a_dx_sph,
									const double a_dx,
									const double a_dy,
									const double a_dz)
	{
		double R_t = (inputs.r_out * c_AU - inputs.r_in * c_AU) / (exp(inputs.C_rad) - 1.0);
		double eta_here[3];
		double eta_ahead[3];
		for (int i = 0; i < DIM; i++)
		{
			double dxd;
			if (i == 0)
				dxd = a_dx;
			if (i == 1)
				dxd = a_dy;
			if (i == 2)
				dxd = a_dz;
			eta_here[i] = a_pt[i] * dxd;
			eta_ahead[i] = a_pt[i] * dxd + dxd;
		}

		double r_here = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_here[0]) - 1.0);
		double theta_here = c_PI * eta_here[1];
		double phi_here = 2.0 * c_PI * eta_here[2];

		double r_ahead = inputs.r_in * c_AU + R_t * (exp(inputs.C_rad * eta_ahead[0]) - 1.0);
		double theta_ahead = c_PI * eta_ahead[1];
		double phi_ahead = 2.0 * c_PI * eta_ahead[2];

		double dr = r_ahead - r_here;
		double dtheta = theta_ahead - theta_here;
		double dphi = phi_ahead - phi_here;

		a_dx_sph(0) = dr;	  // dr
		a_dx_sph(1) = dtheta; // dtheta
		a_dx_sph(2) = dphi;	  // dphi
	}
	PROTO_KERNEL_END(get_delta_sph_coords_calcF, get_delta_sph_coords_calc)

	void get_delta_sph_coords(BoxData<double, DIM> &a_dx_sph,
							  const Box &a_bx,
							  const double a_dx,
							  const double a_dy,
							  const double a_dz)
	{
		forallInPlace_p(get_delta_sph_coords_calc, a_bx, a_dx_sph, a_dx, a_dy, a_dz);
	}
	PROTO_KERNEL_START
	void Correct_V_theta_phi_at_poles_calcF(const Point &a_pt,
											Var<double, NUMCOMPS> &a_U_Sph_ave,
											const double a_dx,
											const double a_dy,
											const double a_dz)
	{
		double E2 = (a_pt[1] + 0.5) * a_dy;
		if (E2 < 0.0 || E2 > 1.0)
		{
			a_U_Sph_ave(2) *= -1.0;
			a_U_Sph_ave(3) *= -1.0;

			a_U_Sph_ave(6) *= -1.0;
			a_U_Sph_ave(7) *= -1.0;
		}
	}
	PROTO_KERNEL_END(Correct_V_theta_phi_at_poles_calcF, Correct_V_theta_phi_at_poles_calc)

	void Correct_V_theta_phi_at_poles(BoxData<double, NUMCOMPS> &a_U_Sph_ave,
									  const double a_dx,
									  const double a_dy,
									  const double a_dz)
	{
		forallInPlace_p(Correct_V_theta_phi_at_poles_calc, a_U_Sph_ave, a_dx, a_dy, a_dz);
	}

	void Spherical_2O_map_filling_func(MHDLevelDataState &a_state)
	{
		for (auto dit : a_state.m_cell_volume)
		{
			Box dbx0 = a_state.m_cell_volume[dit].box();
			MHD_Mapping::get_cell_volume(a_state.m_cell_volume[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz);
			MHD_Mapping::get_face_area(a_state.m_face_area[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz);
			MHD_Mapping::get_delta_sph_coords(a_state.m_dx_sph[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz);
			MHD_Mapping::get_sph_coords_cc(a_state.m_x_sph_cc[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz);
			MHD_Mapping::get_sph_coords_fc(a_state.m_x_sph_fc_1[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz, 0);
			MHD_Mapping::get_sph_coords_fc(a_state.m_x_sph_fc_2[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz, 1);
			MHD_Mapping::get_sph_coords_fc(a_state.m_x_sph_fc_3[dit], dbx0, a_state.m_dx, a_state.m_dy, a_state.m_dz, 2);
		}
	}
}
