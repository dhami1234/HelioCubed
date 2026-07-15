#include "MHD_Gravity.H"

#include "MHD_Constants.H"

namespace MHD_Gravity
{
	PROTO_KERNEL_START
	void AddGravitySourceF(
	        const Point& a_pt,
	        State& a_Rhs,
	        const State& a_W_sph,
	        const V& a_x_sph,
	        double a_inv_velocity_scale_sq)
	{
		const double radius = a_x_sph(0);
		const double theta = a_x_sph(1);
		const double phi = a_x_sph(2);
		const double radial_acceleration =
		        -c_G*c_M0/(radius*radius)*a_inv_velocity_scale_sq;
		const double rho = a_W_sph(iRHO);

		// Momentum is stored in Cartesian components.
		a_Rhs(iMOMX) += rho*radial_acceleration*sin(theta)*cos(phi);
		a_Rhs(iMOMY) += rho*radial_acceleration*sin(theta)*sin(phi);
		a_Rhs(iMOMZ) += rho*radial_acceleration*cos(theta);

		// Gravitational work changes the total MHD energy.
		a_Rhs(iE) += rho*a_W_sph(iVX)*radial_acceleration;
	}
	PROTO_KERNEL_END(AddGravitySourceF, AddGravitySource)

	void add_gravity_source(
	        BoxData<double,NUMCOMPS>& a_Rhs,
	        const BoxData<double,NUMCOMPS>& a_W_sph,
	        const BoxData<double,DIM>& a_x_sph,
	        double a_velocity_scale)
	{
		const double inv_velocity_scale_sq =
		        1.0/(a_velocity_scale*a_velocity_scale);
		forallInPlace_p(
		        AddGravitySource, a_Rhs, a_W_sph, a_x_sph,
		        inv_velocity_scale_sq);
	}
}
