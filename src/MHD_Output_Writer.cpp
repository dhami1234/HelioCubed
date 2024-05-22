#include "Proto.H"
#include "MHD_Output_Writer.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
// #include "Proto_WriteBoxData.H"
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
					const double dt,
					bool CME_data)
	{
		double dx = state.m_dx;
		double dy = state.m_dy;
		double dz = state.m_dz;
		LevelBoxData<double,NUMCOMPS> new_state(state.m_dbl,Point::Ones(NGHOST));
		LevelBoxData<double,NUMCOMPS> new_state2(state.m_dbl,Point::Ones(NGHOST));
		LevelBoxData<double,NUMCOMPS> new_state3(state.m_dbl,Point::Zeros());
		LevelBoxData<double,DIM> phys_coords(state.m_dbl,Point::Ones(NGHOST));
		LevelBoxData<double,NUMCOMPS+DIM> out_data(state.m_dbl,Point::Ones(NGHOST));
		LevelBoxData<double,NUMCOMPS> out_data2(state.m_dbl,Point::Zeros());
		LevelBoxData<double,DIM> x_sph(state.m_dbl,Point::Zeros());

		for (auto dit : new_state){		
				
			state.m_U[dit].copyTo(new_state2[dit]);
			MHDOp::NonDimToDimcalc(new_state2[dit]);
			MHDOp::consToPrimcalc(new_state3[ dit],new_state2[ dit],inputs.gamma);
			// MHDOp::consToPrimcalc(new_state3[ dit],state.m_U[ dit],inputs.gamma);
			if (inputs.output_in_spherical_coords == 1){
				MHD_Mapping::get_sph_coords_cc(x_sph[ dit],x_sph[ dit].box(),dx, dy, dz);
				MHD_Mapping::Cartesian_to_Spherical(new_state[ dit],new_state3[ dit],x_sph[ dit]);
			} else {
				new_state3[ dit].copyTo(new_state[ dit]);
			}

			MHD_Mapping::phys_coords_calc(phys_coords[ dit],state.m_U[ dit].box(),dx,dy,dz);
			MHD_Mapping::out_data_calc(out_data[ dit],phys_coords[ dit],new_state[ dit]);
		}
	
		std::string filename_Data=inputs.Data_file_Prefix+std::to_string(k);
		if (CME_data){
			filename_Data=inputs.Data_file_Prefix+"at_CME_insertion_"+std::to_string(k);
		}
		HDF5Handler h5;
		h5.setTime(time);
		h5.setTimestep(dt);

		if (inputs.output_in_spherical_coords == 1){
			#if TURB == 1
			h5.writeLevel({"X","Y","Z","density","Vr","Vt","Vp", "p","Br","Bt","Bp","Z2","Sigma_c","Lambda"}, 1, out_data, filename_Data);
			#else
			h5.writeLevel({"X","Y","Z","density","Vr","Vt","Vp", "p","Br","Bt","Bp"}, 1, out_data, filename_Data);
			#endif
		} else {
			#if TURB == 1
			h5.writeLevel({"X","Y","Z","density","Vx","Vy","Vz", "p","Bx","By","Bz","Z2","Sigma_c","Lambda"}, 1, out_data, filename_Data);
			#else
			h5.writeLevel({"X","Y","Z","density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, out_data, filename_Data);
			#endif
		}
		
		if(procID()==0) cout << "Written data file after step "<< k << endl;	
	}


	void Write_checkpoint(MHDLevelDataState& state,
					const int k,
					const double time,
					const double dt,
					bool CME_checkpoint)
	{
		double dx = state.m_dx;
		double dy = state.m_dy;
		double dz = state.m_dz;
		LevelBoxData<double,NUMCOMPS> out_data2(state.m_dbl,Point::Zeros());
		std::string filename_Checkpoint=inputs.Checkpoint_file_Prefix+std::to_string(k);
		if (CME_checkpoint){
			filename_Checkpoint=inputs.Checkpoint_file_Prefix+"before_CME_"+std::to_string(k);
		}
		(state.m_U).copyTo(out_data2);
		HDF5Handler h5;
		h5.setTime(time);
		h5.setTimestep(dt);
		
		#if TURB == 1
		h5.writeLevel({"density","Vx","Vy","Vz", "p","Bx","By","Bz", "Z2","Sigma_c","Lambda"}, 1, out_data2, filename_Checkpoint);
		#else
		h5.writeLevel({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, out_data2, filename_Checkpoint);
		#endif
	
		if(procID()==0) cout << "Written checkpoint file after step "<< k << endl;	
	}
}
