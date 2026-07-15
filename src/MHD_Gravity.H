#pragma once

#include "Proto.H"
#include "MHDOp.H"

namespace MHD_Gravity
{
	/**
	 * Add the solar-gravity source to a nondimensional Cartesian conserved RHS.
	 * The primitive state and coordinates are expressed in spherical components.
	 */
	void add_gravity_source(
	        BoxData<double,NUMCOMPS>& a_Rhs,
	        const BoxData<double,NUMCOMPS>& a_W_sph,
	        const BoxData<double,DIM>& a_x_sph,
	        double a_velocity_scale);
}
