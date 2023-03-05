#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include "Proto.H"
#include "MHDLevelDataRK4.H"
// #include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
#include "MHD_Initialize.H"
#include "MHD_EulerStep.H"
#include "MHDOp.H"
#include "MHD_Mapping.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Constants.H"
#include "MHD_CFL.H"
#include "MHD_Probe.H"
#include "MHD_Pre_Time_Step.H"
#include <chrono> // Used by timer
#include "MHDReader.H"
#include "RK4.H"
#include "PolarExchangeCopier.H"
#include "MHD_Set_Boundary_Values.H"

using namespace std;
using namespace Proto;
using namespace MHD_EulerStep;

Parsefrominputs inputs;
int main(int argc, char* argv[])
{
#ifdef PR_MPI
	MPI_Init(&argc,&argv);
#endif
	//have to do this to get a time table
	PR_TIMER_SETFILE("proto.time.table");
	PR_TIME("main");
	int pid = procID();
	//Reading inputs file
	inputs.parsenow(argc,argv);
	#ifdef PR_MPI
		MPI_Barrier(MPI_COMM_WORLD);
	#endif
	int maxLev;	
	// When using mapping, computational domain is always from 0 to 1. The physical grid is mapped from this cube.
	if (inputs.grid_type_global > 1){
		inputs.domsizex = 1.0;
		inputs.domsizey = 1.0;
		inputs.domsizez = 1.0;
	}
	bool takeviscositystep = false;
	if ((inputs.non_linear_visc_apply == 1) || (inputs.linear_visc_apply == 1)) takeviscositystep = true;

	LevelBoxData<double,NUMCOMPS> U[3];  // Size 3 is needed for the convergence rate tests (If indicated in inputs file)
	if (inputs.convTestType != 0) {
		maxLev = 3;
	} else {
		maxLev = 1;
	}
	for (int lev=0; lev<maxLev; lev++)
	{
		// Creating a box for our full domain 
		#if DIM == 2
			Box domain(Point::Zeros(),Point(inputs.domainSizex-1, inputs.domainSizey-1));
		#endif
		#if DIM == 3
			Box domain(Point::Zeros(),Point(inputs.domainSizex-1, inputs.domainSizey-1, inputs.domainSizez-1));
		#endif
		array<bool,DIM> per;
		// All outer boundaries are set to periodic by default
		for(int idir = 0; idir < DIM; idir++){
			per[idir]=true;
	 	}
		// Creating problem domain 
		ProblemDomain pd(domain,per);
		double dx = inputs.domsizex/inputs.domainSizex, dy = inputs.domsizey/inputs.domainSizey, dz = inputs.domsizez/inputs.domainSizez;
		double dt;
		// Following is done for required dt control in convergence tests
		if (inputs.convTestType == 1)
		{
			dt = inputs.CFL*inputs.velocity_scale;
		} else {
			#if DIM == 2
			dt = inputs.CFL*std::min({dx,dy})*inputs.velocity_scale;
			#endif
			#if DIM == 3
			dt = inputs.CFL*std::min({dx,dy,dz})*inputs.velocity_scale;
			#endif
		}
		if (inputs.convTestType == 2)
		{
			dt /= pow(2,lev);
		}
		if (inputs.convTestType == 0) dt = 0.;
		// Create an object state. state.m_U has all the consereved variables (multiplied by Jacobian for mapped grids)
		// All the mapping variables, which are functions of mapping geometry are also included in this class object.
		MHDLevelDataState state(pd,inputs.BoxSize*Point::Ones(),dx, dy, dz, inputs.gamma);
		(state.m_U).setToZero();  

		// This is used to find number of boxes in each processor.
		int count=0;
		for (auto dit : state.m_U)
		{
			count++;
		}
		std::cout << "proc_id: " << pid << ";      num boxes: " << count << std::endl;

		if (inputs.grid_type_global == 2){
				MHD_Mapping::Spherical_2O_map_filling_func(state);
				MHD_Mapping::Spherical_map_filling_func(state);
				MHD_Mapping::Spherical_map_filling_func2(state);
		} else {
			MHD_Mapping::Regular_map_filling_func(state);
		}
		
		double time = 0.;
		double dt_new = 0.;
		double time_seconds;
		MHDReader reader;
		HDF5Handler h5;

		// Read data from h5 BC
		std::vector<BoxData<double, NUMCOMPS, HOST>> BC_data;

		// std::vector<double> dtheta;
		if (inputs.grid_type_global == 2) reader.readData(BC_data, inputs.BC_file);

		// if (inputs.grid_type_global == 2) reader.readGeom(dtheta, inputs.BC_file);
		// if (procID()) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, BC_data[0], "OFT_BCs");
		
		state.m_CME_inserted = false;
		state.m_CME_checkpoint_written = false;

		if (inputs.restartStep == 0){
			if (inputs.grid_type_global == 2 && (inputs.initialize_in_spherical_coords == 1)){
				if (inputs.Spherical_2nd_order == 0) MHD_Initialize::initializeState_Spherical(state);
				if (inputs.Spherical_2nd_order == 1) MHD_Initialize::initializeState_Spherical_2O(state);
			} else {
				MHD_Initialize::initializeState(state);
			}
		} else {
			std::string filename_Checkpoint=inputs.Checkpoint_file_Prefix+std::to_string(inputs.restartStep);
			LevelBoxData<double,NUMCOMPS> readData(state.m_dbl,Point::Zero()); 
			h5.readLevel(readData, filename_Checkpoint);
			for (auto dit : state.m_U){	
				(readData[ dit]).copyTo(state.m_U[ dit]);
			}
			time = h5.time();
			dt = h5.dt();
			time *=inputs.velocity_scale;
			dt *=inputs.velocity_scale;
			double physical_time = MHD_Probe::getPhysTime(time/inputs.velocity_scale);
			if (physical_time > inputs.CME_Enter_Time){
				state.m_CME_inserted = true;
				state.m_CME_checkpoint_written = true;
			}
		}
		
		int start_iter = 0;
		if (inputs.restartStep != 0) {start_iter = inputs.restartStep;}
		if(pid==0) cout << "starting time loop from step " << start_iter << " , maxStep = " << inputs.maxStep << endl;
		bool give_space_in_probe_file = true;
		double probe_cadence = 0;
		double dt_old = dt;
		

		for (int k = start_iter; (k <= inputs.maxStep) && (time/inputs.velocity_scale < inputs.tstop); k++)
		{	
			auto start = chrono::steady_clock::now();
			state.m_divB_calculated = false;
			state.m_divV_calculated = false;
			state.m_Viscosity_calculated = false;
			state.m_min_dt_calculated = false;
			
			for (auto dit : state.m_U){	
				MHDOp::Fix_negative_P(state.m_U[ dit],inputs.gamma); // Current version only for 2nd order spherical	
			}
			if (k!=start_iter){
				if (k!=start_iter+1){
					if (inputs.convTestType == 0){
						dt = inputs.CFL*state.m_min_dt;
						if ((inputs.tstop*inputs.velocity_scale - time) < dt) dt = inputs.tstop*inputs.velocity_scale - time;
					}
				}
				// Below objects need to be created inside the time loop. Otherwise the init in dx keeps on eating memory
				// This is used to take rk4 step
				RK4<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> rk4;
				// This will be used to take Euler step (Primarily used in convergence tests and debugging)
				EulerStep<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> eulerstep;
				// Both Powell divergence cleaning and viscosity implementation need Euler steps at each state update
				EulerStep<MHDLevelDataState, MHDLevelDatadivBOp, MHDLevelDataDX> divBstep;
				EulerStep<MHDLevelDataState, MHDLevelDataViscosityOp, MHDLevelDataDX> viscositystep;

				if (inputs.grid_type_global == 2) MHD_Set_Boundary_Values::interpolate_h5_BC(state, BC_data, time/inputs.velocity_scale);

				double physical_time = MHD_Probe::getPhysTime(time/inputs.velocity_scale);
				double checkpoint_phys_time = inputs.CME_Enter_Time - (inputs.CME_get_checkpoint_before/365.25/24/3600);
				if (!state.m_CME_checkpoint_written && (physical_time >= checkpoint_phys_time)) {
					if (procID() == 0) cout << "Writing Checkpoint " << (inputs.CME_Enter_Time-physical_time)*365.25*24*3600 <<" (s) before CME insertion" << endl;
					MHD_Output_Writer::Write_checkpoint(state, k, time/inputs.velocity_scale, dt/inputs.velocity_scale, true);	
					state.m_CME_checkpoint_written = true;
				}
				MHD_Pre_Time_Step::Insert_CME(state, k, time/inputs.velocity_scale, dt/inputs.velocity_scale);
				if (inputs.convTestType == 1 || inputs.timeIntegratorType == 1) {
					PR_TIME("eulerstep");
					eulerstep.advance(time,dt,state);
				} else {
					if (inputs.timeIntegratorType == 4){
						rk4.advance(time,dt,state);
					}
				}

				if (takeviscositystep) {
					// Take step for artificial viscosity
					viscositystep.advance(time,dt_old,state);
				}

				if (inputs.takedivBstep == 1) {
					// Take step for divB term
					PR_TIME("divBstep");
					divBstep.advance(time,dt,state);
				}

				
				time += dt;
				dt_old = dt;
				time_seconds = time/inputs.velocity_scale;
			}
			
			if (inputs.convTestType == 0)
			{

				int probe_cadence_new = floor(time/inputs.velocity_scale/inputs.probe_cadence);
				if (probe_cadence_new > probe_cadence){
					MHD_Probe::Probe(state,time/inputs.velocity_scale,give_space_in_probe_file);
					give_space_in_probe_file = false;
					probe_cadence = probe_cadence_new;
					if(pid==0) cout << "Probed" << endl;
				}

				if(((inputs.outputInterval > 0) && ((k)%inputs.outputInterval == 0)) || time/inputs.velocity_scale == inputs.tstop || ((inputs.outputInterval > 0) && (k == 0 || k == inputs.restartStep)))
				{	
					if (inputs.sph_inner_BC_hdf5 == 1){
						MHD_Output_Writer::Write_data(state, k, MHD_Probe::getPhysTime(time/inputs.velocity_scale), dt/inputs.velocity_scale, false);
					} else {
						MHD_Output_Writer::Write_data(state, k, time/inputs.velocity_scale, dt/inputs.velocity_scale, false);
					}
								
				}
				if((((inputs.CheckpointInterval > 0) && ((k)%inputs.CheckpointInterval == 0)) || time/inputs.velocity_scale == inputs.tstop || ((inputs.CheckpointInterval > 0) && (k == 0))) && (k!=start_iter || k==0))
				{
					MHD_Output_Writer::Write_checkpoint(state, k, time/inputs.velocity_scale, dt/inputs.velocity_scale, false);	
					std::string filename_to_delete=inputs.Checkpoint_file_Prefix+std::to_string(k-(inputs.CheckpointInterval*inputs.MaxCheckpointFiles))+".hdf5";
					const char* str = filename_to_delete.c_str();
					if (procID() == 0) std::remove(str);
				}
			}
			auto end = chrono::steady_clock::now();
			
			if (inputs.sph_inner_BC_hdf5 == 1){
				double physical_time = MHD_Probe::getPhysTime(time/inputs.velocity_scale);
				if(pid==0) cout <<"nstep = " << k << " Phys Time = " << setprecision(8) << physical_time << " t(s) = " << time/inputs.velocity_scale << " dt(s) = " << dt/inputs.velocity_scale << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms"  << endl;
			} else {
				if(pid==0) cout <<"nstep = " << k << " t(s) = " << time/inputs.velocity_scale << " dt(s) = " << dt/inputs.velocity_scale << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms"  << endl;
			}
		}	
	
		if (inputs.convTestType != 0) {
			//Solution on a single patch
			#if DIM == 2
			U[lev].define(DisjointBoxLayout(pd,Point(inputs.domainSizex, inputs.domainSizey)),Point::Zeros());
			#endif
			#if DIM == 3
			U[lev].define(DisjointBoxLayout(pd,Point(inputs.domainSizex, inputs.domainSizey, inputs.domainSizez)), Point::Zeros());
			#endif
			(state.m_U).copyTo(U[lev]);
			if (inputs.convTestType != 2){
				inputs.domainSizex *= 2;
				inputs.domainSizey *= 2;
				inputs.domainSizez *= 2;
				inputs.BoxSize *= 2; //For debugging: if you want to keep the number of boxes the same
			}
			if (inputs.convTestType == 2 || inputs.convTestType == 3){
				inputs.maxStep *= 2;
			}
		}
	}
	//Until we have a coarsening operation for LevelBoxData, we perform the error calculations on a single patch.
	if(pid==0 && (inputs.convTestType != 0))
	{
		for (int varr = 0; varr < NUMCOMPS; varr++) {
			double ErrMax[2];
			for(int ilev=0; ilev<2; ilev++)
			{
				auto dit_lev=U[ilev].begin();
				auto dit_levp1=U[ilev+1].begin();

				BoxData<double,1> err=slice(U[ilev][*dit_lev],varr);
				err-=Stencil<double>::AvgDown(2)(slice(U[ilev+1][*dit_levp1],varr));
				ErrMax[ilev]=err.absMax();
				std::string filename="Comp_"+std::to_string(varr)+"_err_"+std::to_string(ilev);
				//NOTE: this assumes that the domain length is 1.0, which is assumed throughout this code. May cause errors if this changes.
				double dx=1./(err.box().size(0));
				if (inputs.saveConvTestData){
					HDF5Handler h5;
					h5.writePatch({"err"}, 1, err, filename.c_str());					
				}
				std::cout << "Lev: " << ilev << " , " << ErrMax[ilev] << std::endl;
			}
			double rate = log(abs(ErrMax[0]/ErrMax[1]))/log(2.0);
			std::cout << "order of accuracy = " << rate << std::endl;
		}
	}
	PR_TIMER_REPORT();
#ifdef PR_MPI
	MPI_Finalize();
#endif
}
