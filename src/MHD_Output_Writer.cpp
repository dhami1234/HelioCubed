#include "Proto.H"
#include "MHD_Output_Writer.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
#include "Proto_LevelBoxData.H"
#include "Proto_ProblemDomain.H"
#include "MHD_Mapping.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Constants.H"
#include "MHDLevelDataRK4.H"
extern Parsefrominputs inputs;
/// @brief MHD_Output_Writer namespace
namespace MHD_Output_Writer {

	void Write_data(MHDLevelDataState& state,
					const int k,
					const double time,
					const double dt)
	{
		double dx = state.m_dx;
		double dy = state.m_dy;
		double dz = state.m_dz;
		LevelBoxData<double,NUMCOMPS> new_state(state.m_dbl,Point::Ones(NGHOST));
		LevelBoxData<double,NUMCOMPS> new_state3(state.m_dbl,Point::Zeros());
		LevelBoxData<double,DIM> phys_coords(state.m_dbl,Point::Ones(NGHOST));
		LevelBoxData<double,NUMCOMPS+DIM> out_data(state.m_dbl,Point::Ones(NGHOST));
		LevelBoxData<double,NUMCOMPS> out_data2(state.m_dbl,Point::Zeros());
		LevelBoxData<double,DIM> x_sph(state.m_dbl,Point::Zeros());

		for (auto dit : new_state){		
			if (inputs.grid_type_global == 2){
				if (inputs.Spherical_2nd_order == 0){
					MHD_Mapping::JU_to_W_Sph_ave_calc_func(new_state[ dit], state.m_U[ dit], (state.m_detAA_inv_avg)[ dit], (state.m_A_inv_avg)[ dit], (state.m_r2rdot_avg)[ dit], (state.m_detA_avg)[ dit], (state.m_A_row_mag_avg)[ dit], inputs.gamma, true);
				}
				if (inputs.Spherical_2nd_order == 1){
					MHDOp::consToPrimcalc(new_state3[ dit],state.m_U[ dit],inputs.gamma);
					MHD_Mapping::get_sph_coords_cc(x_sph[ dit],x_sph[ dit].box(),dx, dy, dz);
					MHD_Mapping::Cartesian_to_Spherical(new_state[ dit],new_state3[ dit],x_sph[ dit]);
				}
			} else {
				//W_bar itself is not 4th order W. But it is calculated from 4th order accurate JU for output.
				MHD_Mapping::JU_to_W_bar_calc(new_state[ dit],state.m_U[ dit],(state.m_J)[ dit],inputs.gamma);
			}
			MHD_Mapping::phys_coords_calc(phys_coords[ dit],state.m_U[ dit].box(),dx,dy,dz);
			MHD_Mapping::out_data_calc(out_data[ dit],phys_coords[ dit],new_state[ dit]);
		}
	
		std::string filename_Data=inputs.Data_file_Prefix+std::to_string(k);
		HDF5Handler h5;
		h5.setTime(time);
		h5.setTimestep(dt);
		#if DIM == 2
		h5.writeLevel({"X","Y","density","Vx","Vy", "p","Bx","By"}, 1, out_data, filename_Data);
		#endif
		#if DIM == 3
		h5.writeLevel({"X","Y","Z","density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, out_data, filename_Data);
		#endif
		if(procID()==0) cout << "Written data file after step "<< k << endl;	
	}


	void Write_checkpoint(MHDLevelDataState& state,
					const int k,
					const double time,
					const double dt)
	{
		double dx = state.m_dx;
		double dy = state.m_dy;
		double dz = state.m_dz;
		LevelBoxData<double,NUMCOMPS> out_data2(state.m_dbl,Point::Zeros());
		std::string filename_Checkpoint=inputs.Checkpoint_file_Prefix+std::to_string(k);
		(state.m_U).copyTo(out_data2);
		HDF5Handler h5;
		h5.setTime(time);
		h5.setTimestep(dt);
		#if DIM == 2
		h5.writeLevel({"density","Vx","Vy", "p","Bx","By"}, 1, out_data2, filename_Checkpoint);
		#endif
		#if DIM == 3
		h5.writeLevel({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, out_data2, filename_Checkpoint);
		#endif
		if(procID()==0) cout << "Written checkpoint file after step "<< k << endl;	
	}
}
