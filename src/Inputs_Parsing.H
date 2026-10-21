#pragma once

#include "Proto.H"
#include "InputParser.H"
#include <algorithm> // For std::remove
#include <iostream>  // For std::pout(0), etc.
#include <iomanip>   // For std::setw, etc.
#include <vector>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include "MHD_pout.H"
#include "MHD_Constants.H"
#include "MHD_Time.H"

#ifndef GIT_HASH_RAW
#define GIT_HASH_RAW unknown_hash
#endif

#ifndef GIT_MESSAGE_RAW
#define GIT_MESSAGE_RAW unknown_message
#endif

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

#define GIT_HASH    STRINGIFY(GIT_HASH_RAW)
#define GIT_MESSAGE STRINGIFY(GIT_MESSAGE_RAW)

using namespace Proto;
using namespace std;

/// @brief ParseInputs Singleton class
class ParseInputs {
private:
    int convTestType = 0; // 0: No convergence test, 1: Space Convergence test, 2: Space and Time Convergence test
    // -1: run all convergence levels in one MPI launch (legacy behavior);
    // 0, 1, 2: run only that level; 3: compare previously written levels.
    int convergence_level = -1;
    // A positive value overrides the convergence timestep. This lets separately
    // launched levels use the same timestep sequence as the legacy driver.
    // Keep this as text because the shared InputParser currently converts
    // numeric inputs through single precision before assigning to a double.
    string convergence_dt = "-1";
    // 0: file-driven solar wind; 1: analytic solar wind;
    // 2: radial pulse; 3: radial pulse with constant Cartesian B;
    // 4: non-radial pulse; 5: non-radial pulse with constant Cartesian B;
    // 6: strong spherical shock; 7: strong spherical shock with constant B;
    // 8: off-center hydro blast; 9: off-center blast with constant B = 0.1 G.
    // For 8/9 use gamma=1.4 (see exec/inputs_Blast).
    int init_condition_type = 0;
    int domainSize; // Number of cells in non radial direction
    int thickness; // Number of cells in radial direction
    // Required runtime geometry. Parse with stod because InputArgs converts
    // double inputs through single precision.
    string inner_radius;
    string outer_radius;
    string radial_stretch;
    int boxSize_nonrad; // Size of the box for domain decomposition
    int boxSize_rad; // Size of the box for domain decomposition
    int max_iter; 
    double max_time; // Maximum simulation time in seconds
    double CFL;
    int limiter; // 0: unlimited, 1: O4 limiter (minmod at O2), 2: minmod, 3: MC (O2 only)
    // 0: preserve spherical symmetry (default); 1: interpolate Cartesian U.
    int state_conversion = 0;
    // Optional seam treatments; retain the earlier, faster path by default.
    int BC_conserved_averages = 0; // 1: centered initialization and averaged BCs.
    int seam_pressure_correction = 0; // 1: match angular ghost energy to pressure.
    int spatial_order = 4; // 2: common PLM operator, 4: standard/Berta operator
    int spatial_operator = 1; // 1: Existing Euler operator, 2: Berta WENOZ operator
    int Riemann_solver_type = 1; // 1: Rusanov, 2: Roe 8-wave (standard and Berta)
    double entropy_fix_coeff = 0.3; // FullSphere entropy-fix coefficient; 0 disables it
    double Berta_eta_d = 0.02; // Berta detector threshold
    int Berta_detector_type = 0; // 0: Jameson detector, 1: Berta derivative-ratio detector
    int Berta_fallback_type = 0; // 0: WENO3 fallback, 1: PLM fallback
    double strength; // Strength of perturbation for test problems.
    int write_cadence = 10000000000; 
    int slice_cadence = 10000000000;
    double write_time_cadence = 1000000000; // In seconds
    double slice_time_cadence = 1000000000; // In seconds
    int checkpoint_cadence = 1000000000;
    int P_floor_cadence = 1000000000; // Iteration cadence for P_floor
    string data_file_prefix = "data";
    string checkpoint_file_prefix = "checkpoint";
    std::string slices_raw;
    bool slice_X = false;
    bool slice_Y = false;
    bool slice_Z = false;
    bool slice_R = false;
    void parseSlices() {
        slice_X = slice_Y = slice_Z = slice_R = false;
        std::istringstream iss(slices_raw);
        std::string token;
        while (iss >> token) {
            if (token == "X") slice_X = true;
            else if (token == "Y") slice_Y = true;
            else if (token == "Z") slice_Z = true;
            else if (token == "R") slice_R = true;
        }
    }
    int output_coordinates = 0; // Volume/slice vectors: 0=spherical, 1=Cartesian.
    int scale_data_with_r = 0; // 0: No scaling, 1: Scale with r
    int output_divB = 0; // 0: Do not write div(B), 1: Append div(B) to data output
    int max_checkpoint_files = 1;
    string restart_file;
    double gamma = 1.66667; // Ratio of specific heats
    int temporal_order = 4; // Temporal order of accuracy
    int divergence_cleaning = 1; // 0: No divergence cleaning, 1: Divergence cleaning
    int linear_viscosity = 1; // 0: No linear viscosity, 1: Linear viscosity
    int nonlinear_viscosity = 1; // 0: No nonlinear viscosity, 1: Nonlinear viscosity
    int viscosity_order = 5; // Linear viscosity order: 3 or 5
    int scale_viscosity_with_r = 1; // 0: historical scaling, 1: radius-scaled viscosity
    int boundary_reflux = 1; // 0: don't reflux at block boundaries, 1: reflux
    int radial_refinement = 0; // 0: refine in all directions by 2 for convergence studies;
                           // 1 refine only in radial direction.
    int MBInterp_define;   // old define (0) new define (1)
    // Mapped multiblock AMR. Each ratio list has one value (broadcast) or one per
    // level interface. Both angular coordinates share the nonradial ratio.
    int amr_max_level = 0;
    string amr_ref_ratio_rad = "2";
    string amr_ref_ratio_nonrad = "2";
    string amr_boxSize_rad = "0";
    string amr_boxSize_nonrad = "0";
    int amr_subcycling = 1;
    int amr_verbosity = 1; // 0: coarse steps, 1: level summaries, 2: individual events.
    int amr_regrid_interval = 4;
    double amr_regrid_time_cadence = 0;
    string amr_regrid_time_cadence_raw;
    int amr_buffer = 2;
    // Legacy fallback for density/pressure when their individual options are
    // omitted. Negative thresholds disable a criterion; zero tags any change.
    double amr_threshold = 0.05;
    string amr_refine_density_thresh_raw, amr_refine_pressure_thresh_raw;
    string amr_refine_B_thresh_raw, amr_refine_V_thresh_raw;
    double amr_refine_density_thresh = 0.05;
    double amr_refine_pressure_thresh = 0.05;
    double amr_refine_B_thresh = -1;
    double amr_refine_V_thresh = -1;
    // Finest level a criterion may create. -1 uses the hierarchy maximum;
    // zero disables that criterion without changing its threshold.
    int amr_refine_density_maxlevel = -1;
    int amr_refine_pressure_maxlevel = -1;
    int amr_refine_B_maxlevel = -1;
    int amr_refine_V_maxlevel = -1;
    string amr_refine_regions;
    string BC_file;
    int num_probes = 0;                                
    std::vector<std::string> probe_filenames;           
    std::vector<std::string> probe_coords;  // Accepts HGI coordinates              
    double Probe_cadence = 1000000000; // In seconds
    double BC_start_time = 0.0; // Like 2024.0 for year 2024
    // Epochs must bypass InputArgs' single-precision numeric conversion.
    string BC_start_time_raw = "0";
    CalendarEpoch BC_start_epoch;
    double BC_cadence = 1.0; // In hours
    int BC_frame_rotate = 0; // 0: Frames don't rotate (SWiG), 1: Frames rotate
    int BC_force_single_frame = 0; // 1: Use only one frame for BCs, 0: Use all frames for BCs (if present)
    int Sun_gravity = 0; // 0: No gravity, 1: Include Sun's gravity as source term
    int CME_type = 0; // 0: No CME, 1: CME
    double CME_Enter_Time = 0.0; // Decimal year; 0 inserts at simulation start.
    string CME_Enter_Time_raw = "0";
    double CME_insertion_time = 0.0; // Seconds after BC_start_time.
    double CME_r1;
    double CME_r0;
    double CME_halfangle;
    double CME_pol_flux_control;
    double CME_FRiED_n;
    double CME_lat;
    double CME_lon;
    double CME_tilt;
    double CME_Vandas_R0;
    double CME_apex_speed;
    double CME_energy_control;
    double CME_density;
    double CME_helicity_sign;
    double CME_a1;
    double CME_shift;
    
    // Private Constructor
    ParseInputs() {}

    // Private Copy Constructor and Assignment Operator
    ParseInputs(const ParseInputs&) = delete;
    ParseInputs& operator=(const ParseInputs&) = delete;

public:
    // Static method to get the instance of the class
    static ParseInputs& getInstance() {
        static ParseInputs instance; // Guaranteed to be destroyed and instantiated on first use
        return instance;
    }

    // Static Getters
    static int get_convTestType() {return getInstance().convTestType;}
    static int get_convergence_level() {return getInstance().convergence_level;}
    static double get_convergence_dt() {return std::stod(getInstance().convergence_dt);}
    static int get_init_condition_type() {return getInstance().init_condition_type;}
    static int get_domainSize() {return getInstance().domainSize;}
    static int get_thickness() {return getInstance().thickness;}
    static double get_inner_radius() {return cubedSphereShellGeometry().r0;}
    static double get_outer_radius() {return cubedSphereShellGeometry().r1;}
    static double get_radial_stretch() {return cubedSphereShellGeometry().radialStretch;}
    static int get_boxSize_nonrad() {return getInstance().boxSize_nonrad;}
    static int get_boxSize_rad() {return getInstance().boxSize_rad;}
    static int get_max_iter() {return getInstance().max_iter;}
    static double get_max_time() {return getInstance().max_time;}
    static double get_CFL() {return getInstance().CFL;}
    static int get_limiter() {return getInstance().limiter;}
    static int get_state_conversion() {return getInstance().state_conversion;}
    static int get_BC_conserved_averages() {return getInstance().BC_conserved_averages;}
    static int get_seam_pressure_correction() {return getInstance().seam_pressure_correction;}
    static int get_spatial_operator() {return getInstance().spatial_operator;}
    static int get_Riemann_solver_type() {return getInstance().Riemann_solver_type;}
    static double get_entropy_fix_coeff() {return getInstance().entropy_fix_coeff;}
    static double get_Berta_eta_d() {return getInstance().Berta_eta_d;}
    static int get_Berta_detector_type() {return getInstance().Berta_detector_type;}
    static int get_Berta_fallback_type() {return getInstance().Berta_fallback_type;}
    static double get_strength() {return getInstance().strength;}
    static int get_write_cadence() {return getInstance().write_cadence;}
    static int get_slice_cadence() {return getInstance().slice_cadence;}
    static double get_write_time_cadence() {return getInstance().write_time_cadence;}
    static double get_slice_time_cadence() {return getInstance().slice_time_cadence;}
    static int get_checkpoint_cadence() {return getInstance().checkpoint_cadence;}
    static int get_P_floor_cadence() {return getInstance().P_floor_cadence;}
    static string get_data_file_prefix() {return getInstance().data_file_prefix;}
    static string get_checkpoint_file_prefix() {return getInstance().checkpoint_file_prefix;}
    static bool get_slice_X() { return getInstance().slice_X; }
    static bool get_slice_Y() { return getInstance().slice_Y; }
    static bool get_slice_Z() { return getInstance().slice_Z; }
    static bool get_slice_R() { return getInstance().slice_R; }
    static int get_output_coordinates() {return getInstance().output_coordinates;}
    static int get_scale_data_with_r() {return getInstance().scale_data_with_r;}
    static int get_output_divB() {return getInstance().output_divB;}
    static int get_max_checkpoint_files() {return getInstance().max_checkpoint_files;}
    static string get_restart_file() {return getInstance().restart_file;}
    static double get_gamma() {return getInstance().gamma;}
    static int get_spatial_order() {return getInstance().spatial_order;}
    static int get_amr_max_level() {return getInstance().amr_max_level;}
    static const string& get_amr_ref_ratio_rad() {return getInstance().amr_ref_ratio_rad;}
    static const string& get_amr_ref_ratio_nonrad() {return getInstance().amr_ref_ratio_nonrad;}
    static const string& get_amr_boxSize_rad() {return getInstance().amr_boxSize_rad;}
    static const string& get_amr_boxSize_nonrad() {return getInstance().amr_boxSize_nonrad;}
    static int get_amr_subcycling() {return getInstance().amr_subcycling;}
    static int get_amr_verbosity() {return getInstance().amr_verbosity;}
    static int get_amr_regrid_interval() {return getInstance().amr_regrid_interval;}
    static double get_amr_regrid_time_cadence() {return getInstance().amr_regrid_time_cadence;}
    static int get_amr_buffer() {return getInstance().amr_buffer;}
    static double get_amr_threshold() {return getInstance().amr_threshold;}
    static double get_amr_refine_density_thresh() {return getInstance().amr_refine_density_thresh;}
    static double get_amr_refine_pressure_thresh() {return getInstance().amr_refine_pressure_thresh;}
    static double get_amr_refine_B_thresh() {return getInstance().amr_refine_B_thresh;}
    static double get_amr_refine_V_thresh() {return getInstance().amr_refine_V_thresh;}
    static int get_amr_refine_density_maxlevel()
    {return getInstance().amr_refine_density_maxlevel < 0 ? get_amr_max_level() : getInstance().amr_refine_density_maxlevel;}
    static int get_amr_refine_pressure_maxlevel()
    {return getInstance().amr_refine_pressure_maxlevel < 0 ? get_amr_max_level() : getInstance().amr_refine_pressure_maxlevel;}
    static int get_amr_refine_B_maxlevel()
    {return getInstance().amr_refine_B_maxlevel < 0 ? get_amr_max_level() : getInstance().amr_refine_B_maxlevel;}
    static int get_amr_refine_V_maxlevel()
    {return getInstance().amr_refine_V_maxlevel < 0 ? get_amr_max_level() : getInstance().amr_refine_V_maxlevel;}
    static bool uses_amr_tagging()
    {
        return (get_amr_refine_density_thresh() >= 0 && get_amr_refine_density_maxlevel() > 0)
            || (get_amr_refine_pressure_thresh() >= 0 && get_amr_refine_pressure_maxlevel() > 0)
            || (get_amr_refine_B_thresh() >= 0 && get_amr_refine_B_maxlevel() > 0)
            || (get_amr_refine_V_thresh() >= 0 && get_amr_refine_V_maxlevel() > 0);
    }
    static const string& get_amr_refine_regions() {return getInstance().amr_refine_regions;}
    static int get_temporal_order() {return getInstance().temporal_order;}
    static int get_divergence_cleaning() {return getInstance().divergence_cleaning;}
    static bool uses_second_order_limiter()
    {return get_spatial_order() == 2 && get_limiter() != 0;}
    // The slope limiter supplies dissipation for limited O2. Apply the same switches
    // to fluxes, timestep estimation, and ghost-width selection.
    static int get_linear_viscosity()
    {return uses_second_order_limiter() ? 0 : getInstance().linear_viscosity;}
    static int get_nonlinear_viscosity()
    {return uses_second_order_limiter() ? 0 : getInstance().nonlinear_viscosity;}
    static int get_viscosity_order() {return getInstance().viscosity_order;}
    static int get_scale_viscosity_with_r() {return getInstance().scale_viscosity_with_r;}
    static int get_boundary_reflux() {return getInstance().boundary_reflux;}
    static int get_radial_refinement() {return getInstance().radial_refinement;}
    static int get_MBInterp_define() {return getInstance().MBInterp_define;}
    static string get_BC_file() {return getInstance().BC_file;}
    static int get_num_probes() { return getInstance().num_probes; }                               
    static const std::vector<std::string>& get_probe_filenames() { return getInstance().probe_filenames; } 
    static const std::vector<std::string>& get_probe_coords()     { return getInstance().probe_coords; }     
    static double get_Probe_cadence() {return getInstance().Probe_cadence;}
    static double get_BC_start_time() {return getInstance().BC_start_time;}
    static const CalendarEpoch& get_BC_start_epoch() {return getInstance().BC_start_epoch;}
    static double get_BC_cadence() {return getInstance().BC_cadence;}
    static int get_BC_frame_rotate() {return getInstance().BC_frame_rotate;}
    static int get_BC_force_single_frame() {return getInstance().BC_force_single_frame;}
    static int get_Sun_gravity() {return getInstance().Sun_gravity;}
    static int get_CME_type() {return getInstance().CME_type;}
    static double get_CME_Enter_Time() {return getInstance().CME_Enter_Time;}
    static double get_CME_insertion_time() {return getInstance().CME_insertion_time;}
    static double get_CME_r1() {return getInstance().CME_r1;}
    static double get_CME_r0() {return getInstance().CME_r0;}
    static double get_CME_halfangle() {return getInstance().CME_halfangle;}
    static double get_CME_pol_flux_control() {return getInstance().CME_pol_flux_control;}
    static double get_CME_FRiED_n() {return getInstance().CME_FRiED_n;}
    static double get_CME_lat() {return getInstance().CME_lat;}
    static double get_CME_lon() {return getInstance().CME_lon;}
    static double get_CME_tilt() {return getInstance().CME_tilt;}
    static double get_CME_Vandas_R0() {return getInstance().CME_Vandas_R0;}
    static double get_CME_apex_speed() {return getInstance().CME_apex_speed;}
    static double get_CME_energy_control() {return getInstance().CME_energy_control;}
    static double get_CME_density() {return getInstance().CME_density;}
    static double get_CME_helicity_sign() {return getInstance().CME_helicity_sign;}
    static double get_CME_a1() {return getInstance().CME_a1;}
    static double get_CME_shift() {return getInstance().CME_shift;}

    // Parses now with command-line arguments
    void parsenow(int argc, char* argv[]) {
        InputArgs args;
        args.add("convTestType", convTestType);
        args.add("convergence_level", convergence_level);
        args.add("convergence_dt", convergence_dt);
        args.add("init_condition_type", init_condition_type);
        args.add("domainSize", domainSize);
        args.add("thickness", thickness);
        args.add("inner_radius", inner_radius);
        args.add("outer_radius", outer_radius);
        args.add("radial_stretch", radial_stretch);
        args.add("boxSize_nonrad", boxSize_nonrad);
        args.add("boxSize_rad", boxSize_rad);
        args.add("max_iter", max_iter);
        args.add("max_time", max_time);
        args.add("CFL", CFL);
        args.add("limiter", limiter);
        args.add("state_conversion", state_conversion);
        args.add("BC_conserved_averages", BC_conserved_averages);
        args.add("seam_pressure_correction", seam_pressure_correction);
        args.add("spatial_order", spatial_order);
        args.add("spatial_operator", spatial_operator);
        args.add("Riemann_solver_type", Riemann_solver_type);
        args.add("entropy_fix_coeff", entropy_fix_coeff);
        args.add("Berta_eta_d", Berta_eta_d);
        args.add("Berta_detector_type", Berta_detector_type);
        args.add("Berta_fallback_type", Berta_fallback_type);
        args.add("strength", strength);
        args.add("write_cadence", write_cadence);
        args.add("slice_cadence", slice_cadence);
        args.add("write_time_cadence", write_time_cadence);
        args.add("slice_time_cadence", slice_time_cadence);
        args.add("checkpoint_cadence", checkpoint_cadence);
        args.add("P_floor_cadence", P_floor_cadence);
        args.add("data_file_prefix", data_file_prefix);
        args.add("checkpoint_file_prefix", checkpoint_file_prefix);
        args.add("slices", slices_raw);
        args.add("output_coordinates", output_coordinates);
        args.add("scale_data_with_r", scale_data_with_r);
        args.add("output_divB", output_divB);
        args.add("max_checkpoint_files", max_checkpoint_files);
        args.add("restart_file", restart_file);
        args.add("gamma", gamma);
        args.add("temporal_order", temporal_order);
        args.add("divergence_cleaning", divergence_cleaning);
        args.add("linear_viscosity", linear_viscosity);
        args.add("nonlinear_viscosity", nonlinear_viscosity);
        args.add("viscosity_order", viscosity_order);
        args.add("scale_viscosity_with_r", scale_viscosity_with_r);
        args.add("boundary_reflux", boundary_reflux);
        args.add("radial_refinement", radial_refinement);
        args.add("MBInterp_define", MBInterp_define);
        args.add("amr_max_level", amr_max_level);
        args.add("amr_ref_ratio_rad", amr_ref_ratio_rad);
        args.add("amr_ref_ratio_nonrad", amr_ref_ratio_nonrad);
        args.add("amr_boxSize_rad", amr_boxSize_rad);
        args.add("amr_boxSize_nonrad", amr_boxSize_nonrad);
        args.add("amr_subcycling", amr_subcycling);
        args.add("amr_verbosity", amr_verbosity);
        args.add("amr_regrid_interval", amr_regrid_interval);
        args.add("amr_regrid_time_cadence", amr_regrid_time_cadence_raw);
        args.add("amr_buffer", amr_buffer);
        args.add("amr_threshold", amr_threshold);
        args.add("amr_refine_density_thresh", amr_refine_density_thresh_raw);
        args.add("amr_refine_pressure_thresh", amr_refine_pressure_thresh_raw);
        args.add("amr_refine_B_thresh", amr_refine_B_thresh_raw);
        args.add("amr_refine_V_thresh", amr_refine_V_thresh_raw);
        args.add("amr_refine_density_maxlevel", amr_refine_density_maxlevel);
        args.add("amr_refine_pressure_maxlevel", amr_refine_pressure_maxlevel);
        args.add("amr_refine_B_maxlevel", amr_refine_B_maxlevel);
        args.add("amr_refine_V_maxlevel", amr_refine_V_maxlevel);
        args.add("amr_refine_regions", amr_refine_regions);
        args.add("BC_file", BC_file);
        args.add("num_probes", num_probes); 
        constexpr int MAX_PROBES = 32; // or more if needed
        // ensure storage exists before parse:
        probe_filenames.assign(MAX_PROBES, "");
        probe_coords.assign(MAX_PROBES, "");
        for (int i = 0; i < MAX_PROBES; ++i) {
            const std::string base = "probe" + std::to_string(i) + ".";
            args.add(base + "filename", probe_filenames[(size_t)i]);
            args.add(base + "coords",   probe_coords[(size_t)i]);
        }
        args.add("Probe_cadence", Probe_cadence);
        args.add("BC_start_time", BC_start_time_raw);
        args.add("BC_cadence", BC_cadence);
        args.add("BC_frame_rotate", BC_frame_rotate);
        args.add("BC_force_single_frame", BC_force_single_frame);
        args.add("Sun_gravity", Sun_gravity);
        args.add("CME_Enter_Time", CME_Enter_Time_raw);
        args.add("CME_type", CME_type);
        args.add("CME_r1", CME_r1);
        args.add("CME_r0", CME_r0);
        args.add("CME_halfangle", CME_halfangle);
        args.add("CME_pol_flux_control", CME_pol_flux_control);
        args.add("CME_FRiED_n", CME_FRiED_n);
        args.add("CME_lat", CME_lat);
        args.add("CME_lon", CME_lon);
        args.add("CME_tilt", CME_tilt);
        args.add("CME_Vandas_R0", CME_Vandas_R0);
        args.add("CME_apex_speed", CME_apex_speed);
        args.add("CME_energy_control", CME_energy_control);
        args.add("CME_density", CME_density);
        args.add("CME_helicity_sign", CME_helicity_sign);
        args.add("CME_a1", CME_a1);
        args.add("CME_shift", CME_shift);
        args.parse(argc, argv);
        auto requireGeometry = [](const string& raw, const char* name) {
            if (raw.empty()) {
                throw std::invalid_argument(
                    "HelioCubed requires -inner_radius, -outer_radius, and -radial_stretch "
                    "in the inputs file; -" + string(name) + " is missing or commented out.");
            }
        };
        requireGeometry(inner_radius, "inner_radius");
        requireGeometry(outer_radius, "outer_radius");
        requireGeometry(radial_stretch, "radial_stretch");
        if (spatial_order != 2 && spatial_order != 4) {
            throw std::invalid_argument("spatial_order must be 2 or 4.");
        }
        if (amr_max_level < 0 || amr_regrid_interval < 0
            || amr_buffer < 0 || !std::isfinite(amr_threshold)) {
            throw std::invalid_argument("AMR requires nonnegative max level, regrid interval, and buffer, and a finite threshold.");
        }
        if (amr_subcycling != 0 && amr_subcycling != 1) {
            throw std::invalid_argument("amr_subcycling must be 0 or 1.");
        }
        if (amr_verbosity < 0 || amr_verbosity > 2) {
            throw std::invalid_argument("amr_verbosity must be 0, 1, or 2.");
        }
        if (limiter < 0 || limiter > 3) {
            throw std::invalid_argument("limiter must be 0 (unlimited), 1 (O4 limiter; minmod at O2), 2 (minmod), or 3 (MC at O2).");
        }
        if (limiter == 3 && spatial_order != 2) {
            throw std::invalid_argument("limiter=3 (MC) requires spatial_order=2.");
        }
        if (temporal_order != 1 && temporal_order != 2 && temporal_order != 3 && temporal_order != 4) {
            throw std::invalid_argument("temporal_order must be 1, 2, 3, or 4.");
        }
        if (state_conversion != 0 && state_conversion != 1) {
            throw std::invalid_argument("state_conversion must be 0 (symmetry preserving) or 1 (Cartesian).");
        }
        if (output_coordinates != 0 && output_coordinates != 1) {
            throw std::invalid_argument("output_coordinates must be 0 (spherical) or 1 (Cartesian).");
        }
        if (BC_conserved_averages != 0 && BC_conserved_averages != 1) {
            throw std::invalid_argument("BC_conserved_averages must be 0 or 1.");
        }
        if (seam_pressure_correction != 0 && seam_pressure_correction != 1) {
            throw std::invalid_argument("seam_pressure_correction must be 0 or 1.");
        }
        if ((spatial_order == 2 || state_conversion == 1) && CUBED_SPHERE_SHELL_RADIAL_COORD != 0) {
            throw std::invalid_argument("Cartesian state_conversion requires radial coordinate 0.");
        }
        if (viscosity_order != 3 && viscosity_order != 5) {
            throw std::invalid_argument("viscosity_order must be 3 or 5.");
        }
        auto geometry = cubedSphereShellGeometry();
        auto parsePreciseValue = [](const string& raw, const char* name,
                                     double defaultValue) {
            if (raw.empty()) return defaultValue;
            try {
                size_t end = 0;
                double value = std::stod(raw, &end);
                size_t tail = raw.find_first_not_of(" \t\r\n", end);
                if (tail != string::npos && raw[tail] != '#') {
                    throw std::invalid_argument("trailing characters");
                }
                return value;
            } catch (const std::exception&) {
                throw std::invalid_argument(string(name) + " must be a number.");
            }
        };
        auto refinementThreshold = [&](const string& raw, const char* name, double fallback) {
            const double value = parsePreciseValue(raw, name, fallback);
            if (!std::isfinite(value))
                throw std::invalid_argument(string(name) + " must be finite.");
            return value;
        };
        amr_regrid_time_cadence = parsePreciseValue(amr_regrid_time_cadence_raw,
            "amr_regrid_time_cadence", 0);
        if (!std::isfinite(amr_regrid_time_cadence) || amr_regrid_time_cadence < 0)
            throw std::invalid_argument("amr_regrid_time_cadence must be finite and nonnegative (0 disables it).");
        amr_refine_density_thresh = refinementThreshold(amr_refine_density_thresh_raw,
            "amr_refine_density_thresh", amr_threshold);
        amr_refine_pressure_thresh = refinementThreshold(amr_refine_pressure_thresh_raw,
            "amr_refine_pressure_thresh", amr_threshold);
        amr_refine_B_thresh = refinementThreshold(amr_refine_B_thresh_raw, "amr_refine_B_thresh", -1);
        amr_refine_V_thresh = refinementThreshold(amr_refine_V_thresh_raw, "amr_refine_V_thresh", -1);
        if (amr_refine_density_maxlevel < -1 || amr_refine_pressure_maxlevel < -1
            || amr_refine_B_maxlevel < -1 || amr_refine_V_maxlevel < -1)
            throw std::invalid_argument("AMR variable maxlevels must be -1 (hierarchy maximum) or nonnegative.");
        BC_start_epoch = parseCalendarEpoch(BC_start_time_raw, "BC_start_time");
        CalendarEpoch cmeEpoch = parseCalendarEpoch(CME_Enter_Time_raw, "CME_Enter_Time");
        BC_start_time = BC_start_epoch.decimalYear();
        CME_Enter_Time = cmeEpoch.decimalYear();
        if (CME_type != 0 && CME_type != 1) {
            throw std::invalid_argument("CME_type must be 0 or 1.");
        }
        if (CME_type == 1 && CME_Enter_Time > 0.0 && BC_start_time == 0.0) {
            throw std::invalid_argument(
                "A scheduled CME requires BC_start_time to define the simulation's calendar epoch.");
        }
        CME_insertion_time = CME_type == 1 && cmeEpoch.year != 0
            ? secondsBetweenEpochs(BC_start_epoch, cmeEpoch) : 0.0;
        geometry.r0 = parsePreciseValue(inner_radius, "inner_radius", 0.0);
        geometry.r1 = parsePreciseValue(outer_radius, "outer_radius", 0.0);
        geometry.radialStretch = parsePreciseValue(radial_stretch, "radial_stretch", 0.0);
        if (!std::isfinite(geometry.r0) || !std::isfinite(geometry.r1)
            || geometry.r0 <= 0.0 || geometry.r1 <= geometry.r0) {
            throw std::invalid_argument(
                "Shell radii must be finite and satisfy 0 < inner_radius < outer_radius (cm).");
        }
        if (!std::isfinite(geometry.radialStretch)
            || !std::isfinite(std::expm1(geometry.radialStretch))) {
            throw std::invalid_argument("radial_stretch must be finite without exponential overflow.");
        }
        // Configure Proto before any maps, interpolation operators, or metrics
        // are built, so every consumer uses the same runtime geometry.
        cubedSphereShellGeometry() = geometry;
        parseSlices();
        // use only the first num_probes entries
        if (num_probes < 0 || num_probes > MAX_PROBES)
            throw std::invalid_argument("num_probes must be between 0 and 32.");
        probe_filenames.resize((size_t)std::max(0, num_probes));
        probe_coords.resize((size_t)std::max(0, num_probes));

        pout(0) << setw(75) << setfill('-') << "-" << endl;
        pout(0) << "Last git pull info:" << endl;
        pout(0) << setw(75) << setfill('-') << "-" << endl;
        pout(0) << "Git hash: " << GIT_HASH << endl; 
        pout(0) << "Commit message: " << GIT_MESSAGE << endl;
        pout(0) << setw(75) << setfill('-') << "-" << endl;
        pout(0) << "Running with the following Cubed Sphere geometry from inputs" << endl;
        pout(0) << setw(75) << setfill('-') << "-" << endl;
        pout(0) << "Inner Radius: " << get_inner_radius()/c_SR << " SR" << endl;
        pout(0) << "Outer Radius: " << get_outer_radius()/c_SR << " SR" << endl;
        pout(0) << "Radial Grid Stretching Parameter (C): " << get_radial_stretch() << endl;
        pout(0) << setw(75) << setfill('-') << "-" << endl;
        pout(0) << "Running with the following inputs (these can be changed in the inputs file)" << endl;
        pout(0) << setw(75) << setfill('-') << "-" << endl;
        pout(0) << "convTestType: " << convTestType << endl;
        pout(0) << "convergence_level: " << convergence_level << endl;
        pout(0) << "convergence_dt: " << convergence_dt << endl;
        pout(0) << "init_condition_type: " << init_condition_type << endl;
        pout(0) << "domainSize: " << domainSize << endl;
        pout(0) << "thickness: " << thickness << endl;
        pout(0) << "boxSize_nonrad: " << boxSize_nonrad << endl;
        pout(0) << "boxSize_rad: " << boxSize_rad << endl;
        pout(0) << "max_iter: " << max_iter << endl;
        pout(0) << "max_time: " << max_time << endl;
        pout(0) << "CFL: " << CFL << endl;
        pout(0) << "limiter: " << limiter << endl;
        pout(0) << "state_conversion: " << state_conversion << endl;
        pout(0) << "BC_conserved_averages: " << BC_conserved_averages << endl;
        pout(0) << "seam_pressure_correction: " << seam_pressure_correction << endl;
        pout(0) << "spatial_operator: " << spatial_operator << endl;
        pout(0) << "Riemann_solver_type: " << Riemann_solver_type << endl;
        pout(0) << "entropy_fix_coeff: " << entropy_fix_coeff << endl;
        pout(0) << "Berta_eta_d: " << Berta_eta_d << endl;
        pout(0) << "Berta_detector_type: " << Berta_detector_type << endl;
        pout(0) << "Berta_fallback_type: " << Berta_fallback_type << endl;
        pout(0) << "strength: " << strength << endl;
        pout(0) << "write_cadence: " << write_cadence << endl;
        pout(0) << "slice_cadence: " << slice_cadence << endl;
        pout(0) << "write_time_cadence: " << write_time_cadence << endl;
        pout(0) << "slice_time_cadence: " << slice_time_cadence << endl;
        pout(0) << "checkpoint_cadence: " << checkpoint_cadence << endl;
        pout(0) << "P_floor_cadence: " << P_floor_cadence << endl;
        pout(0) << "data_file_prefix: " << data_file_prefix << endl;
        pout(0) << "checkpoint_file_prefix: " << checkpoint_file_prefix << endl;
        pout(0) << "Slices: " << slices_raw << endl;
        pout(0) << "output_coordinates: " << output_coordinates
                << (output_coordinates == 0 ? " (spherical)" : " (Cartesian)") << endl;
        pout(0) << "scale_data_with_r: " << scale_data_with_r << endl;
        pout(0) << "output_divB: " << output_divB << endl;
        pout(0) << "max_checkpoint_files: " << max_checkpoint_files << endl;
        pout(0) << "restart_file: " << restart_file << endl;
        pout(0) << "gamma: " << gamma << endl;
        pout(0) << "spatial_order: " << spatial_order << endl;
        if (amr_max_level) {
            pout(0) << "AMR: max level " << amr_max_level << ", radial ratios " << amr_ref_ratio_rad
                    << ", angular ratios " << amr_ref_ratio_nonrad << (temporal_order==4 ? ", RK4 with " : temporal_order==3 ? ", RK3 with " : ", SSPRK2 with ")
                    << (amr_subcycling ? "CFL-based time subcycling" : "a common timestep")
                    << "; regrid interval "
                    << amr_regrid_interval << " base steps, time cadence " << amr_regrid_time_cadence
                    << "(s), buffer " << amr_buffer << endl;
            pout(0) << "AMR fine patch sizes: radial " << amr_boxSize_rad
                    << ", angular " << amr_boxSize_nonrad << " (0=automatic)" << endl;
            auto criterion = [&](const char* name, double threshold, int maxLevel) {
                pout(0) << "AMR relative-gradient criterion " << name << ": ";
                if (threshold < 0 || maxLevel == 0) pout(0) << "disabled";
                else pout(0) << "threshold " << threshold << ", finest level " << std::min(maxLevel, amr_max_level);
                pout(0) << endl;
            };
            criterion("density", amr_refine_density_thresh, get_amr_refine_density_maxlevel());
            criterion("pressure", amr_refine_pressure_thresh, get_amr_refine_pressure_maxlevel());
            criterion("|B|", amr_refine_B_thresh, get_amr_refine_B_maxlevel());
            criterion("|v|", amr_refine_V_thresh, get_amr_refine_V_maxlevel());
        }
        if (spatial_order == 2) {
            const char* name = limiter == 0 ? "unlimited" : limiter == 3 ? "MC (monotonized central)" : "minmod";
            pout(0) << "Spatial reconstruction: Cartesian PLM with " << name
                    << "; state_conversion, spatial_operator, Berta settings, and seam_pressure_correction are unused." << endl;
        }
        pout(0) << "temporal_order: " << temporal_order << endl;
        pout(0) << "divergence_cleaning: " << divergence_cleaning << endl;
        pout(0) << "linear_viscosity: " << get_linear_viscosity() << endl;
        pout(0) << "nonlinear_viscosity: " << get_nonlinear_viscosity() << endl;
        if (uses_second_order_limiter() && (linear_viscosity || nonlinear_viscosity))
            pout(0) << "Artificial viscosity inputs are ignored for spatial_order=2 with minmod or MC; the slope limiter supplies the dissipation." << endl;
        pout(0) << "viscosity_order: " << viscosity_order << endl;
        pout(0) << "scale_viscosity_with_r: " << scale_viscosity_with_r << endl;
        pout(0) << "boundary_reflux: " << boundary_reflux << endl;
        pout(0) << "radial_refinement: " << radial_refinement << endl;
        pout(0) << "MBInterp_define: " << MBInterp_define << endl;
        pout(0) << "BC_file: " << BC_file << endl;
        pout(0) << "Probe_cadence: " << Probe_cadence << endl;
        pout(0) << "num_probes: " << num_probes << endl; 
        for (int i = 0; i < num_probes; ++i) {        
            pout(0) << "  probe" << i << ".filename: " << probe_filenames[static_cast<size_t>(i)] << endl;
            pout(0) << "  probe" << i << ".coords:   " << probe_coords[static_cast<size_t>(i)]   << endl;
        }
        pout(0) << "Probe_cadence: " << Probe_cadence << endl;
        pout(0) << "BC_start_time: " << BC_start_time_raw << endl;
        pout(0) << "BC_cadence: " << BC_cadence << endl;
        pout(0) << "BC_frame_rotate: " << BC_frame_rotate << endl;
        pout(0) << "BC_force_single_frame: " << BC_force_single_frame << endl;
        pout(0) << "Sun_gravity: " << Sun_gravity << endl;
        pout(0) << "CME_type: " << CME_type << endl;
        pout(0) << "CME_Enter_Time: " << CME_Enter_Time_raw << endl;
        pout(0) << "CME_r1: " << CME_r1 << endl;
        pout(0) << "CME_r0: " << CME_r0 << endl;
        pout(0) << "CME_halfangle: " << CME_halfangle << endl;
        pout(0) << "CME_pol_flux_control: " << CME_pol_flux_control << endl;
        pout(0) << "CME_FRiED_n: " << CME_FRiED_n << endl;
        pout(0) << "CME_lat: " << CME_lat << endl;
        pout(0) << "CME_lon: " << CME_lon << endl;
        pout(0) << "CME_tilt: " << CME_tilt << endl;
        pout(0) << "CME_Vandas_R0: " << CME_Vandas_R0 << endl;
        pout(0) << "CME_apex_speed: " << CME_apex_speed << endl;
        pout(0) << "CME_energy_control: " << CME_energy_control << endl;
        pout(0) << "CME_density: " << CME_density << endl;
        pout(0) << "CME_helicity_sign: " << CME_helicity_sign << endl;
        pout(0) << "CME_a1: " << CME_a1 << endl;
        pout(0) << "CME_shift: " << CME_shift << endl;
        pout(0) << setw(75) << setfill('-') << "-" << endl;
                
    }

    // Utility function to remove spaces from a string
    string removeSpaces(string str) {
        str.erase(remove(str.begin(), str.end(), ' '), str.end());
        return str;
    }

    // Additional getters and setters as needed
};
