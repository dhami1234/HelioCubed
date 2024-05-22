#pragma once

#include "Proto.H"
#include "MHDOp.H"
#include "MHDLevelDataRK4.H"
//TODO: These variables need to go into a common file or be defined in a more accessible way
//constexpr unsigned int NUMCOMPS=2*DIM+2;
typedef Var<double,NUMCOMPS> State;


namespace MHD_Mapping {
	void eta_calc(BoxData<double,DIM>& a_eta,
	              const Box& a_bx,
	              const double a_dx,
	              const double a_dy,
	              const double a_dz);

	void etaCorner_calc(BoxData<double,DIM>& a_eta,
	              const Box& a_bx,
	              const double a_dx,
	              const double a_dy,
	              const double a_dz);

	void etaFace_calc(BoxData<double,DIM>& a_eta,
	                  const Box& a_bx,
	                  const double a_dx,
	                  const double a_dy,
	                  const double a_dz,
	                  int a_d);

	void eta_to_x_calc(BoxData<double,DIM>& a_x,
	                   const BoxData<double,DIM>& a_eta,
					   const Box& dbx0);
			   

	void phys_coords_calc(BoxData<double,DIM>& a_x,
	                      const Box& dbx1,
	                      const double a_dx,
	                      const double a_dy,
	                      const double a_dz);	

	void phys_coords_face_calc(BoxData<double,DIM>& a_x,
	                      const Box& dbx1,
	                      const double a_dx,
	                      const double a_dy,
	                      const double a_dz,				  		  
	                      const int a_d);		  			  

	void out_data_calc(BoxData<double,NUMCOMPS+DIM>& a_out_data,
	                   const BoxData<double,DIM>& a_phys_coords,
	                   const BoxData<double,NUMCOMPS>& a_W);

	void Cartesian_to_Spherical(BoxData<double,NUMCOMPS>& a_W_Sph,
	                    		const BoxData<double,NUMCOMPS>& a_W_Cart,
	                    		const BoxData<double,DIM>& a_x_Sph);

	void Spherical_to_Cartesian(BoxData<double,NUMCOMPS>& a_W_Cart,
	                    		const BoxData<double,NUMCOMPS>& a_W_Sph,
	                    		const BoxData<double,DIM>& a_x_Sph);

	void get_sph_coords_fc(BoxData<double,DIM>& a_x_sph,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz,
	                    int a_d);														

	void get_sph_coords_cc(BoxData<double,DIM>& a_x_sph,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz);

	void get_cell_volume(BoxData<double,1>& a_V,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz);

	void get_face_area(BoxData<double,DIM>& a_A,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz);	

	void get_delta_sph_coords(BoxData<double,DIM>& a_dx_sph,
	                    const Box& a_bx,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz);	
	
	void Correct_V_theta_phi_at_poles(BoxData<double,NUMCOMPS>& a_U_Sph_ave,
	                             const double a_dx,
	                             const double a_dy,
	                             const double a_dz);				
								 			   
	void Spherical_2O_map_filling_func(MHDLevelDataState& a_state);
}
