#include "Proto.H"
#include "MHDLevelDataRK4.H"
#include "MHD_Initialize.H"
#include "MHD_EulerStep.H"
#include "MHDOp.H"
#include "MHD_Mapping.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Probe.H"
#include "MHD_Pre_Time_Step.H"
#include <chrono> // Used by timer
#include "MHDReader.H"
#include "RK4.H"
#include "MHD_Set_Boundary_Values.H"

using namespace std;
using namespace Proto;
using namespace MHD_EulerStep;

Parsefrominputs inputs;



PROTO_KERNEL_START
void fixBCnumsF(Var<double, NUMCOMPS> &a_out_data,
				Var<double, NUMCOMPS-3> &a_BC)
{
	for (int i = 0; i < NUMCOMPS-3; i++)
	{
			a_out_data(i) = a_BC(i); 
	}
	for (int i = NUMCOMPS-3; i < NUMCOMPS; i++)
	{
			a_out_data(i) = 0.; 
	}
}
PROTO_KERNEL_END(fixBCnumsF, fixBCnums)


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
	
	// Creating a box for our full domain 
	Box domain(Point::Zeros(),Point(inputs.domainSizex-1, inputs.domainSizey-1, inputs.domainSizez-1));
	
	array<bool,DIM> per;
	// All outer boundaries are set to periodic by default
	for(int idir = 0; idir < DIM; idir++){
		per[idir]=true;
	}
	// Creating problem domain 
	ProblemDomain pd(domain,per);
	double dx = 1.0/inputs.domainSizex, dy = 1.0/inputs.domainSizey, dz = 1.0/inputs.domainSizez;
	double dt = 0.;
	
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


	MHD_Mapping::Spherical_2O_map_filling_func(state);

	
	
	double time = 0.;
	double dt_new = 0.;
	double time_seconds;
	MHDReader reader;
	HDF5Handler h5;
	std::vector<BoxData<double, NUMCOMPS, HOST>> BC_data;
	// Read data from h5 BC
	#if TURB == 1
		std::vector<BoxData<double, NUMCOMPS-3, HOST>> BC_data_temp;
		if (inputs.sph_inner_BC_hdf5 == 1) reader.readData(BC_data_temp, inputs.BC_file);
		int siz_h5 = BC_data_temp.size();
		BC_data.resize(siz_h5);
		for (int i =0; i<siz_h5; i++){
			BC_data[i] = forall<double, NUMCOMPS>(fixBCnums, BC_data_temp[i]);
		}
	#else
		if (inputs.sph_inner_BC_hdf5 == 1) reader.readData(BC_data, inputs.BC_file);
	#endif

	
	// std::vector<double> dtheta;
	// if (inputs.grid_type_global == 2) reader.readGeom(dtheta, inputs.BC_file);
	// if (procID()) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, BC_data[0], "OFT_BCs");
	
	state.m_CME_inserted = false;
	state.m_CME_checkpoint_written = false;

	if (inputs.restartStep == 0){

		MHD_Initialize::initializeState_Spherical_2O(state);

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
		}
		double checkpoint_phys_time = inputs.CME_Enter_Time - (inputs.CME_get_checkpoint_before/365.25/24/3600);
		if (!state.m_CME_checkpoint_written && (physical_time >= checkpoint_phys_time)) {
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
		state.m_min_dt_calculated = false;
		
		for (auto dit : state.m_U){	
			MHDOp::Fix_negative_P(state.m_U[ dit],inputs.gamma); // Current version only for 2nd order spherical	
		}
		if (k!=start_iter){
			if (k!=start_iter+1){
				
				dt = inputs.CFL*state.m_min_dt;
				if ((inputs.tstop*inputs.velocity_scale - time) < dt) dt = inputs.tstop*inputs.velocity_scale - time;
				
			}
			// Below objects need to be created inside the time loop. Otherwise the init in dx keeps on eating memory
			// This is used to take rk4 step
			RK4<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> rk4;
			// This will be used to take Euler step (Primarily used in convergence tests and debugging)
			EulerStep<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> eulerstep;
			// Both Powell divergence cleaning and viscosity implementation need Euler steps at each state update
			EulerStep<MHDLevelDataState, MHDLevelDatadivBOp, MHDLevelDataDX> divBstep;

			if (inputs.sph_inner_BC_hdf5 == 1) MHD_Set_Boundary_Values::interpolate_h5_BC(state, BC_data, time/inputs.velocity_scale);

			
			MHD_Pre_Time_Step::Insert_CME(state, k, time/inputs.velocity_scale, dt/inputs.velocity_scale);
			if (inputs.timeIntegratorType == 1) {
				eulerstep.advance(time,dt,state);
			} else {
				if (inputs.timeIntegratorType == 4){
					rk4.advance(time,dt,state);
				}
			}

			if (inputs.takedivBstep == 1) {
				// Take step for divB term
				PR_TIME("divBstep");
				divBstep.advance(time,dt,state);
			}

			
			time += dt;
			dt_old = dt;
			time_seconds = time/inputs.velocity_scale;
			state.m_time = time;
			double physical_time = MHD_Probe::getPhysTime(time/inputs.velocity_scale);
			double checkpoint_phys_time = inputs.CME_Enter_Time - (inputs.CME_get_checkpoint_before/365.25/24/3600);
			if (!state.m_CME_checkpoint_written && (physical_time >= checkpoint_phys_time)) {
				if (procID() == 0) cout << "Writing Checkpoint " << (inputs.CME_Enter_Time-physical_time)*365.25*24*3600 <<" (s) before CME insertion" << endl;
				MHD_Output_Writer::Write_checkpoint(state, k, time/inputs.velocity_scale, dt/inputs.velocity_scale, true);	
				state.m_CME_checkpoint_written = true;
			}
		}
		
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
		
		auto end = chrono::steady_clock::now();
		
		if (inputs.sph_inner_BC_hdf5 == 1){
			double physical_time = MHD_Probe::getPhysTime(time/inputs.velocity_scale);
			if(pid==0) cout <<"nstep = " << k << " Phys Time = " << setprecision(8) << physical_time << " t(s) = " << time/inputs.velocity_scale << " dt(s) = " << dt/inputs.velocity_scale << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms"  << endl;
		} else {
			if(pid==0) cout <<"nstep = " << k << " t(s) = " << time/inputs.velocity_scale << " dt(s) = " << dt/inputs.velocity_scale << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms"  << endl;
		}
	}	
	
	PR_TIMER_REPORT();
#ifdef PR_MPI
	MPI_Finalize();
#endif
}
