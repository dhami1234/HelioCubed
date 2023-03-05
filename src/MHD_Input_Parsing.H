#pragma once

#include "Proto.H"
#include "InputParser.H"

using namespace Proto;
/// @brief Parsefrominputs class
class Parsefrominputs {
  public:
    double tstop;
	double CFL;
	double domsizex;
	double domsizey;
	double domsizez;
	double gamma;
	int domainSizex;
	int domainSizey;
	int domainSizez;
	int maxStep;
	int restartStep;
	int outputInterval;
	int CheckpointInterval;
	int MaxCheckpointFiles;
    int BoxSize;
    int grid_type_global; // 0: 2D-Rectangular/3D-Rectangular;  1: 2D-Wavy/3D-Not implemented;  2: 2D-Polar/3D-Spherical
    int LowBoundType;  // 0 for periodic, 1 for Dirichlet, 2 for open. This is for dir==0 only
    int HighBoundType;  // 0 for periodic, 1 for Dirichlet, 2 for open. This is for dir==0 only
    int sph_inner_BC_hdf5; // 1 is h5 file is used to provide inner BCs
    int limiter_apply;
    int slope_flattening_apply;
    int non_linear_visc_apply;
    int linear_visc_apply;
    int takedivBstep;
    int timeIntegratorType; // 1: Euler Step    4: RK4 
    int init_condition_type;
    /*
	   0. constant solution
	   1. 2D current sheet problem
	   2. flow from 1 side
	   3. 2D Orszag Tang problem
	   4. Alfven wave problem
	   5. Acoustic pulse problem
	   6. Acoustic pulse problem in polar grid
	   7. Cylindrical pulse problem in polar grid
	   8. Shifted Acoustic pulse problem
	   9. Euler problem
	   10. Radially out flow in polar grid
	   11. Cylindrical pulse problem in cartesian grid
	   12. Acoustic pulse problem in 3D cartesian grid
	   13. 3D MHD blast wave
	   14. 2D MHD blast wave
	   15. Acoustic pulse problem with Bx
	   16. Shell acoustic pulse problem in spherical grid
	   17. Sphere Shifted Acoustic pulse problem in Spherical grid

	   19. Velocity pulse problem in polar grid
	   21. Radially out flow in spherical grid
	 */
    int Riemann_solver_type;  //1: Rusanov Solver, 2: Roe Solver
    double entropy_fix_coeff;
    int initialize_in_spherical_coords;
    int output_in_spherical_coords;
    double C1_fix; // A constant in wavy grid definition.
    double r_in;
    double r_out;
    double C_rad; // A constant in exponential dr in spherical grid.
    int convTestType; // 0 for no convergence test, 1 for space convergence, 2 for time convergence, 3 for space and time convergence
	int saveConvTestData;
	int Spherical_2nd_order; // 1 if we want to turn off Phil's mapping and use spherical grid method of MS-FLUKSS
	double BC_start_time; // In decimal year
	double BC_cadence; // In hours
	double probe_cadence; // In seconds
	string BC_file;
	string Data_file_Prefix;
	string Checkpoint_file_Prefix;
	string Probe_trajectory_file;
	string Probe_data_file;
	double CME_Enter_Time;
	double CME_get_checkpoint_before; // Get a checkpoint file this much time in seconds before CME insertion.
	double CME_r1, CME_r0, CME_halfangle, CME_pol_flux_control, CME_FRiED_n, CME_lat, CME_lon, CME_tilt, CME_Vandas_R0, CME_apex_speed, CME_energy_control, CME_density, CME_helicity_sign, CME_a1, CME_shift;
	//MHD equations can be written in dimensionless variables. We use scaling for this.
	double density_scale, velocity_scale;
  public:

	string removeSpaces(string str)
	{
		str.erase(remove(str.begin(), str.end(), ' '), str.end());
		return str;
	}
    
	void parsenow(int argc, char* argv[]){
		InputArgs args;
		args.add("tstop",  tstop);
		args.add("CFL",  CFL);
		args.add("domsizex",  domsizex);
		args.add("domsizey",  domsizey);
		args.add("domsizez",  domsizez);
		args.add("gamma",  gamma);
		args.add("domainSizex",  domainSizex);
		args.add("domainSizey",  domainSizey);
		args.add("domainSizez",  domainSizez);
		args.add("maxStep",  maxStep);
		args.add("restartStep",  restartStep);
		args.add("outputInterval",  outputInterval);
		args.add("CheckpointInterval",  CheckpointInterval);
		args.add("MaxCheckpointFiles",  MaxCheckpointFiles);
		args.add("BoxSize",  BoxSize);
		args.add("grid_type_global",  grid_type_global);
		args.add("LowBoundType",  LowBoundType);
		args.add("HighBoundType",  HighBoundType);
		args.add("sph_inner_BC_hdf5",  sph_inner_BC_hdf5);
		args.add("limiter_apply",  limiter_apply);
		args.add("slope_flattening_apply",  slope_flattening_apply);
		args.add("non_linear_visc_apply",  non_linear_visc_apply);
		args.add("linear_visc_apply",  linear_visc_apply);
		args.add("takedivBstep",  takedivBstep);
		args.add("timeIntegratorType",  timeIntegratorType);
		args.add("init_condition_type",  init_condition_type);
		args.add("Riemann_solver_type",  Riemann_solver_type);
		args.add("entropy_fix_coeff",  entropy_fix_coeff);
		args.add("initialize_in_spherical_coords",  initialize_in_spherical_coords);
		args.add("output_in_spherical_coords",  output_in_spherical_coords);
		args.add("C1_fix",  C1_fix);
		args.add("r_in",  r_in);
		args.add("r_out",  r_out);
		args.add("C_rad",  C_rad);
		args.add("convTestType",  convTestType);
		args.add("saveConvTestData",  saveConvTestData);
		args.add("Spherical_2nd_order",  Spherical_2nd_order);
		args.add("BC_start_time",  BC_start_time);
		args.add("BC_cadence",  BC_cadence);
		args.add("probe_cadence",  probe_cadence);
		args.add("BC_file",  BC_file);
		args.add("Data_file_Prefix",  Data_file_Prefix);
		args.add("Checkpoint_file_Prefix",  Checkpoint_file_Prefix);
		args.add("Probe_trajectory_file",  Probe_trajectory_file);
		args.add("Probe_data_file",  Probe_data_file);
		args.add("CME_Enter_Time",  CME_Enter_Time);
		args.add("CME_get_checkpoint_before",  CME_get_checkpoint_before);
		args.add("CME_r1",  CME_r1);
		args.add("CME_r0",  CME_r0);
		args.add("CME_halfangle",  CME_halfangle);
		args.add("CME_pol_flux_control",  CME_pol_flux_control);
		args.add("CME_FRiED_n",  CME_FRiED_n);
		args.add("CME_lat",  CME_lat);
		args.add("CME_lon",  CME_lon);
		args.add("CME_tilt",  CME_tilt);
		args.add("CME_Vandas_R0",  CME_Vandas_R0);
		args.add("CME_apex_speed",  CME_apex_speed);
		args.add("CME_energy_control",  CME_energy_control);
		args.add("CME_density",  CME_density);
		args.add("CME_helicity_sign",  CME_helicity_sign);
		args.add("CME_a1",  CME_a1);
		args.add("CME_shift",  CME_shift);
		args.add("density_scale",  density_scale);
		args.add("velocity_scale",  velocity_scale);
		args.parse(argc,argv);
		// args.print(); //It prints, but not in the order I want.

		int pid = procID();

		if (pid == 0){
			cout << "Running with the following inputs (these can be changed in the file 'inputs')" << endl;
			cout << setw(75) << setfill('-') << "-" << endl;
			cout << "tstop: " << tstop << endl;
			cout << "CFL: " << CFL << endl;
			cout << "domsizex: " << domsizex << endl;
			cout << "domsizey: " << domsizey << endl;
			cout << "domsizez: " << domsizez << endl;
			cout << "gamma: " << gamma << endl;
			cout << "domainSizex: "  << domainSizex << endl;
			cout << "domainSizey: "  << domainSizey << endl;
			cout << "domainSizez: "  << domainSizez << endl;
			cout << "maxStep: "    << maxStep << endl;
			cout << "restartStep: "    << restartStep << endl;
			cout << "outputInterval: "    << outputInterval << endl;	
			cout << "CheckpointInterval: "    << CheckpointInterval << endl;	
			cout << "MaxCheckpointFiles: "    << MaxCheckpointFiles << endl;	
			cout << "BoxSize: "    << BoxSize << endl;
			cout << "grid_type_global: "    << grid_type_global << endl;
			cout << "LowBoundType: "    << LowBoundType << endl;
			cout << "HighBoundType: "    << HighBoundType << endl;
			cout << "sph_inner_BC_hdf5: "    << sph_inner_BC_hdf5 << endl;
			cout << "limiter_apply: "    << limiter_apply << endl;
			cout << "slope_flattening_apply: "    << slope_flattening_apply << endl;
			cout << "non_linear_visc_apply: "    << non_linear_visc_apply << endl;
			cout << "linear_visc_apply: "    << linear_visc_apply << endl;
			cout << "takedivBstep: "    << takedivBstep << endl;
			cout << "timeIntegratorType: "    << takedivBstep << endl;
			cout << "init_condition_type: "    << init_condition_type << endl;
			cout << "Riemann_solver_type: "    << Riemann_solver_type << endl;
			cout << "entropy_fix_coeff: "    << entropy_fix_coeff << endl;
			cout << "initialize_in_spherical_coords: "    << initialize_in_spherical_coords << endl;
			cout << "output_in_spherical_coords: "    << output_in_spherical_coords << endl;
			cout << "C1_fix: "    << C1_fix << endl;
			cout << "r_in: "    << r_in << endl;
			cout << "r_out: "    << r_out << endl;
			cout << "C_rad: "    << C_rad << endl;
			cout << "convTestType: "    << convTestType << endl;
			cout << "saveConvTestData: "    << saveConvTestData << endl;
			cout << "Spherical_2nd_order: "    << Spherical_2nd_order << endl;
			cout << "BC_start_time: "    << BC_start_time << endl;
			cout << "BC_cadence: "    << BC_cadence << endl;
			cout << "probe_cadence: "    << probe_cadence << endl;
			cout << "BC_file: "    << BC_file << endl;
			cout << "Data_file_Prefix: "    << Data_file_Prefix << endl;
			cout << "Checkpoint_file_Prefix: "    << Checkpoint_file_Prefix << endl;
			cout << "Probe_trajectory_file: "    << Probe_trajectory_file << endl;
			cout << "Probe_data_file: "    << Probe_data_file << endl;
			cout << "CME_Enter_Time: "    << CME_Enter_Time << endl;
			cout << "CME_get_checkpoint_before: "    << CME_get_checkpoint_before << endl;
			cout << "CME_r1: "    << CME_r1 << endl;
			cout << "CME_r0: "    << CME_r0 << endl;
			cout << "CME_halfangle: "    << CME_halfangle << endl;
			cout << "CME_pol_flux_control: "    << CME_pol_flux_control << endl;
			cout << "CME_FRiED_n: "    << CME_FRiED_n << endl;
			cout << "CME_lat: "    << CME_lat << endl;
			cout << "CME_lon: "    << CME_lon << endl;
			cout << "CME_tilt: "    << CME_tilt << endl;
			cout << "CME_Vandas_R0: "    << CME_Vandas_R0 << endl;
			cout << "CME_apex_speed: "    << CME_apex_speed << endl;
			cout << "CME_energy_control: "    << CME_energy_control << endl;
			cout << "CME_density: "    << CME_density << endl;
			cout << "CME_helicity_sign: "    << CME_helicity_sign << endl;
			cout << "CME_a1: "    << CME_a1 << endl;
			cout << "CME_shift: "    << CME_shift << endl;
			cout << "density_scale: "    << density_scale << endl;
			cout << "velocity_scale: "    << velocity_scale << endl;
			cout << setw(75) << setfill('-') << "-" << endl;
		}

	}
    //Use the following format if encapsulation is needed (if variables needed to be made private) 
    //Getter
	int getBoxSize() {
      return BoxSize;
    }
	// Setter
    void setBoxSize(int s) {
      BoxSize = s;
    }

};