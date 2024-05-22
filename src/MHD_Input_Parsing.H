#pragma once

#include "Proto.H"
#include "InputParser.H"

using namespace Proto;
/// @brief Parsefrominputs class
class Parsefrominputs {
  public:
    double tstop;
	double CFL;
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
    int sph_inner_BC_hdf5; // 1 is h5 file is used to provide inner BCs
    int limiter_apply;
    int takedivBstep;
    int timeIntegratorType; // 1: Euler Step    4: RK4 
    int Riemann_solver_type;  //1: Rusanov Solver, 2: Roe Solver
    double entropy_fix_coeff;
    int initialize_in_spherical_coords;
    int output_in_spherical_coords;
    double r_in;
    double r_out;
    double C_rad; // A constant in exponential dr in spherical grid.
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
	double Sun_Z2, Sun_SigmaC, Sun_Lambda, SigmaD, Scaling_time; // Turbulence parameters for the Sun. Scaling time is the time it takes for Z2 to reach the given value in the inner boundary.
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
		args.add("sph_inner_BC_hdf5",  sph_inner_BC_hdf5);
		args.add("limiter_apply",  limiter_apply);
		args.add("takedivBstep",  takedivBstep);
		args.add("timeIntegratorType",  timeIntegratorType);
		args.add("Riemann_solver_type",  Riemann_solver_type);
		args.add("entropy_fix_coeff",  entropy_fix_coeff);
		args.add("initialize_in_spherical_coords",  initialize_in_spherical_coords);
		args.add("output_in_spherical_coords",  output_in_spherical_coords);
		args.add("r_in",  r_in);
		args.add("r_out",  r_out);
		args.add("C_rad",  C_rad);
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
		args.add("Sun_Z2",  Sun_Z2);
		args.add("Sun_SigmaC",  Sun_SigmaC);
		args.add("Sun_Lambda",  Sun_Lambda);
		args.add("SigmaD",  SigmaD);
		args.add("Scaling_time",  Scaling_time);
		args.parse(argc,argv);
		// args.print(); //It prints, but not in the order I want.

		int pid = procID();

		if (pid == 0){
			cout << "Running with the following inputs (these can be changed in the file 'inputs')" << endl;
			cout << setw(75) << setfill('-') << "-" << endl;
			cout << "tstop: " << tstop << endl;
			cout << "CFL: " << CFL << endl;
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
			cout << "sph_inner_BC_hdf5: "    << sph_inner_BC_hdf5 << endl;
			cout << "limiter_apply: "    << limiter_apply << endl;
			cout << "takedivBstep: "    << takedivBstep << endl;
			cout << "timeIntegratorType: "    << takedivBstep << endl;
			cout << "Riemann_solver_type: "    << Riemann_solver_type << endl;
			cout << "entropy_fix_coeff: "    << entropy_fix_coeff << endl;
			cout << "initialize_in_spherical_coords: "    << initialize_in_spherical_coords << endl;
			cout << "output_in_spherical_coords: "    << output_in_spherical_coords << endl;
			cout << "r_in: "    << r_in << endl;
			cout << "r_out: "    << r_out << endl;
			cout << "C_rad: "    << C_rad << endl;
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
			cout << "Sun_Z2: "    << Sun_Z2 << endl;
			cout << "Sun_SigmaC: "    << Sun_SigmaC << endl;
			cout << "Sun_Lambda: "    << Sun_Lambda << endl;
			cout << "SigmaD: "    << SigmaD << endl;
			cout << "Scaling_time: "    << Scaling_time << endl;
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