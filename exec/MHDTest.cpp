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
#include "MHD_AMR_Euler.H"
#include "AMRRK4.H"
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
#include "MHD_BoxOp.H"
#include "MHD_divB.H"
#include <chrono> // Used by timer
#include "MHDReader.H"
#include "RK4.H"
#include "PolarExchangeCopier.H"
#include "MHD_Set_Boundary_Values.H"

using namespace std;
using namespace Proto;
using namespace MHD_EulerStep;
bool CME_inserted;
Parsefrominputs inputs;
int main(int argc, char *argv[])
{
#ifdef PR_MPI
	MPI_Init(&argc, &argv);
#endif
	// have to do this to get a time table
	PR_TIMER_SETFILE("proto.time.table");
	PR_TIME("main");

	typedef BoxOp_MHD<double> OP;
	typedef BoxOp_divB<double> OP_divB;

	int pid = procID();
	// Reading inputs file
	inputs.parsenow(argc, argv);
#ifdef PR_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	int maxLev;
	// When using mapping, computational domain is always from 0 to 1. The physical grid is mapped from this cube.
	if (inputs.grid_type_global > 1)
	{
		inputs.domsizex = 1.0;
		inputs.domsizey = 1.0;
		inputs.domsizez = 1.0;
	}
	bool takeviscositystep = false;
	if ((inputs.non_linear_visc_apply == 1) || (inputs.linear_visc_apply == 1))
		takeviscositystep = true;

	LevelBoxData<double, NUMCOMPS> Uconv[3]; // Size 3 is needed for the convergence rate tests (If indicated in inputs file)
	if (inputs.convTestType != 0)
	{
		maxLev = 3;
	}
	else
	{
		maxLev = 1;
	}
	for (int lev = 0; lev < maxLev; lev++)
	{
// Creating a box for our full domain
#if DIM == 2
		Box domain(Point::Zeros(), Point(inputs.domainSizex - 1, inputs.domainSizey - 1));
#endif
#if DIM == 3
		Box domain(Point::Zeros(), Point(inputs.domainSizex - 1, inputs.domainSizey - 1, inputs.domainSizez - 1));
#endif
		array<bool, DIM> per;
		// All outer boundaries are set to periodic by default
		for (int idir = 0; idir < DIM; idir++)
		{
			per[idir] = true;
		}
		// Creating problem domain
		ProblemDomain pd(domain, per);

		Point boxSizeVect = Point::Ones(inputs.BoxSize);
		DisjointBoxLayout layout(pd, boxSizeVect);

		double dx = inputs.domsizex / inputs.domainSizex, dy = inputs.domsizey / inputs.domainSizey, dz = inputs.domsizez / inputs.domainSizez;

		Array<double, DIM> dx_arr;
		dx_arr[0] = dx;
#if DIM == 2
		dx_arr[1] = dy;
#endif
#if DIM == 3
		dx_arr[1] = dy;
		dx_arr[2] = dz;
#endif

		double dt;

		// cout << "Creating grid" << endl;
		std::vector<Point> refRatios(inputs.numLevels - 1, Point::Ones(inputs.refRatio));
		AMRGrid grid(layout, refRatios, inputs.numLevels);
		double dxLevel = dx;
		double dyLevel = dy;
		double dzLevel = dz;
		Point bufferSize = Point::Ones(inputs.regridBufferSize);

		if (inputs.numLevels > 1)
		{
			for (int lvl = 0; lvl < inputs.numLevels - 1; lvl++)
			{
				LevelBoxData<double, NUMCOMPS> initData(grid[lvl], Point::Zeros());
				MHD_Initialize::initializeState(initData, dxLevel, dyLevel, dzLevel);
				LevelTagData tags(grid[lvl], bufferSize);
				for (auto iter : grid[lvl])
				{
					OP::generateTags(tags[iter], initData[iter]);
				}
				AMRGrid::buffer(tags, bufferSize);
				grid.regrid(tags, lvl);
				dxLevel /= inputs.refRatio;
				dyLevel /= inputs.refRatio;
				dzLevel /= inputs.refRatio;
			}

			for (int lvl = 2; lvl < inputs.numLevels; lvl++)
			{
				grid.enforceNesting2(lvl);
			}
		}
		// cout << "Grid created" << endl;
		AMRData<double, NUMCOMPS> U(grid, OP::ghost());
		double time = 0.;

		double time_seconds;
		MHDReader reader;
		HDF5Handler h5;

		// Read data from h5 BC
		std::vector<BoxData<double, NUMCOMPS, HOST>> BC_data;

		std::vector<double> dtheta;
		if (inputs.grid_type_global == 2)
			reader.readData(BC_data, inputs.BC_file);

		if (inputs.grid_type_global == 2)
			reader.readGeom(dtheta, inputs.BC_file);
		// if (procID()==0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, BC_data[0], "OFT_BCs");

		CME_inserted = false;
		bool CME_checkpoint_written = false;

		if (inputs.restartStep == 0)
		{
			if (inputs.grid_type_global == 2 && (inputs.initialize_in_spherical_coords == 1))
			{
				// if (inputs.Spherical_2nd_order == 0) MHD_Initialize::initializeState_Spherical(state);
				// if (inputs.Spherical_2nd_order == 1) MHD_Initialize::initializeState_Spherical_2O(state);
			}
			else
			{
				MHD_Initialize::initializeState(U, dx, dy, dz);
			}
		}
		else
		{
			std::string filename_Checkpoint = inputs.Checkpoint_file_Prefix + std::to_string(inputs.restartStep);
			LevelBoxData<double, NUMCOMPS> readData(layout, Point::Zero());
			h5.readLevel(readData, filename_Checkpoint);
			for (auto dit : U[0])
			{
				(readData[dit]).copyTo(U[0][dit]);
			}
			time = h5.time();
			dt = h5.dt();
			time *= inputs.velocity_scale;
			dt *= inputs.velocity_scale;
			double physical_time = MHD_Probe::getPhysTime(time / inputs.velocity_scale);
			if (physical_time > inputs.CME_Enter_Time)
			{
				CME_inserted = true;
			}
			double checkpoint_phys_time = inputs.CME_Enter_Time - (inputs.CME_get_checkpoint_before / 365.25 / 24 / 3600);
			if (!CME_checkpoint_written && (physical_time >= checkpoint_phys_time))
			{
				CME_checkpoint_written = true;
			}
		}
		U.averageDown();
		// cout << "State initialized" << endl;
		// AMRRK4<OP, double, NUMCOMPS> amrrk4(U, dx_arr);
		// AMREuler<OP, double, NUMCOMPS> amreuler(U, dx_arr);
		// AMREuler<OP_divB, double, NUMCOMPS> amrdivB(U, dx_arr);

		AMRRK4<BoxOp_MHD, double> amrrk4(U, dx_arr);
		AMREuler<BoxOp_MHD, double> amreuler(U, dx_arr);
		AMREuler<BoxOp_divB, double> amrdivB(U, dx_arr);

		// cout << "Starting time loop" << endl;

		// This is used to find number of boxes in each processor.
		int count = 0;
		for (int AMRlev = 0; AMRlev < inputs.numLevels; AMRlev++)
			for (auto dit : U[AMRlev])
			{
				count++;
			}
		std::cout << "proc_id: " << pid << ";      num boxes: " << count << std::endl;

		// if (inputs.grid_type_global == 2){
		// 		MHD_Mapping::Spherical_2O_map_filling_func(state);
		// 		MHD_Mapping::Spherical_map_filling_func(state);
		// 		MHD_Mapping::Spherical_map_filling_func2(state);
		// } else {
		// 	MHD_Mapping::Regular_map_filling_func(state);
		// }

		int start_iter = 0;
		if (inputs.restartStep != 0)
		{
			start_iter = inputs.restartStep;
		}
		if (pid == 0)
			cout << "starting time loop from step " << start_iter << " , maxStep = " << inputs.maxStep << endl;
		bool give_space_in_probe_file = true;
		double probe_cadence = 0;
		double dt_old = dt;

		for (int k = start_iter; (k <= inputs.maxStep) && (time / inputs.velocity_scale < inputs.tstop); k++)
		{
			MHD_CFL::dt_calc_func(dt, time, U, dx, dy, dz, lev);
			if (k == start_iter)
				dt = 0.0;

			auto start = chrono::steady_clock::now();
			// state.m_divB_calculated = false;
			// state.m_divV_calculated = false;
			for (int AMRlev = 0; AMRlev < inputs.numLevels; AMRlev++)
				for (auto dit : U[AMRlev])
				{
					MHDOp::Fix_negative_P(U[AMRlev][dit], inputs.gamma); // Current version only for 2nd order spherical
				}

			// 	if (inputs.grid_type_global == 2) MHD_Set_Boundary_Values::interpolate_h5_BC(state, BC_data, time/inputs.velocity_scale);

			MHD_Pre_Time_Step::Insert_CME(U, grid, dx, dy, dz, k, time / inputs.velocity_scale, dt / inputs.velocity_scale);
			if (inputs.convTestType == 1 || inputs.timeIntegratorType == 1)
			{
				PR_TIME("eulerstep");
				amreuler.advance(dt);
			}
			else
			{
				if (inputs.timeIntegratorType == 4)
				{
					amrrk4.advance(dt);
				}
			}

			// 	if (takeviscositystep) {
			// 		// Take step for artificial viscosity
			// 		viscositystep.advance(time,dt_old,state);
			// 	}

			if (inputs.takedivBstep == 1)
			{
				// Take step for divB term
				PR_TIME("divBstep");
				amrdivB.advance(dt);
			}

			time += dt;
			dt_old = dt;
			time_seconds = time / inputs.velocity_scale;

			double physical_time = MHD_Probe::getPhysTime(time / inputs.velocity_scale);
			double checkpoint_phys_time = inputs.CME_Enter_Time - (inputs.CME_get_checkpoint_before / 365.25 / 24 / 3600);
			if (!CME_checkpoint_written && (physical_time >= checkpoint_phys_time))
			{
				if (procID() == 0)
					cout << "Writing Checkpoint " << (inputs.CME_Enter_Time - physical_time) * 365.25 * 24 * 3600 << " (s) before CME insertion" << endl;
				MHD_Output_Writer::Write_checkpoint(U, k, time / inputs.velocity_scale, dt / inputs.velocity_scale, true);
				CME_checkpoint_written = true;
			}

			if (inputs.convTestType == 0)
			{

				int probe_cadence_new = floor(time / inputs.velocity_scale / inputs.probe_cadence);
				if (probe_cadence_new > probe_cadence)
				{
					// MHD_Probe::Probe(state,time/inputs.velocity_scale,give_space_in_probe_file);
					give_space_in_probe_file = false;
					probe_cadence = probe_cadence_new;
					if (pid == 0)
						cout << "Probed" << endl;
				}

				if (((inputs.outputInterval > 0) && ((k) % inputs.outputInterval == 0)) || time / inputs.velocity_scale == inputs.tstop || ((inputs.outputInterval > 0) && (k == 0 || k == inputs.restartStep)))
				{
					if (inputs.sph_inner_BC_hdf5 == 1)
					{
						MHD_Output_Writer::Write_data(U, grid, dx, dy, dz, k, MHD_Probe::getPhysTime(time / inputs.velocity_scale), dt / inputs.velocity_scale, false);
					}
					else
					{
						MHD_Output_Writer::Write_data(U, grid, dx, dy, dz, k, time / inputs.velocity_scale, dt / inputs.velocity_scale, false);
					}
				}
				if ((((inputs.CheckpointInterval > 0) && ((k) % inputs.CheckpointInterval == 0)) || time / inputs.velocity_scale == inputs.tstop || ((inputs.CheckpointInterval > 0) && (k == 0))) && (k != start_iter || k == 0))
				{
					MHD_Output_Writer::Write_checkpoint(U, k, time / inputs.velocity_scale, dt / inputs.velocity_scale, false);
					std::string filename_to_delete = inputs.Checkpoint_file_Prefix + std::to_string(k - (inputs.CheckpointInterval * inputs.MaxCheckpointFiles)) + ".hdf5";
					const char *str = filename_to_delete.c_str();
					if (procID() == 0)
						std::remove(str);
				}
			}
			auto end = chrono::steady_clock::now();

			if (inputs.sph_inner_BC_hdf5 == 1)
			{
				double physical_time = MHD_Probe::getPhysTime(time / inputs.velocity_scale);
				if (pid == 0)
					cout << "nstep = " << k << " Phys Time = " << setprecision(8) << physical_time << " t(s) = " << time / inputs.velocity_scale << " dt(s) = " << dt / inputs.velocity_scale << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
			}
			else
			{
				if (pid == 0)
					cout << "nstep = " << k << " t(s) = " << time / inputs.velocity_scale << " dt(s) = " << dt / inputs.velocity_scale << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
			}
		}

		if (inputs.convTestType != 0)
		{
// Solution on a single patch
#if DIM == 2
			Uconv[lev].define(DisjointBoxLayout(pd, Point(inputs.domainSizex, inputs.domainSizey)), Point::Zeros());
#endif
#if DIM == 3
			Uconv[lev].define(DisjointBoxLayout(pd, Point(inputs.domainSizex, inputs.domainSizey, inputs.domainSizez)), Point::Zeros());
#endif
			(U[0]).copyTo(Uconv[lev]);
			if (inputs.convTestType != 2)
			{
				inputs.domainSizex *= 2;
				inputs.domainSizey *= 2;
				inputs.domainSizez *= 2;
				inputs.BoxSize *= 2; // For debugging: if you want to keep the number of boxes the same
			}
			if (inputs.convTestType == 2 || inputs.convTestType == 3)
			{
				inputs.maxStep *= 2;
			}
		}
		// cout << "lev = " << lev << endl;
	}
	// Until we have a coarsening operation for LevelBoxData, we perform the error calculations on a single patch.
	if (pid == 0 && (inputs.convTestType != 0))
	{
		for (int varr = 0; varr < NUMCOMPS; varr++)
		{
			double ErrMax[2];
			for (int ilev = 0; ilev < 2; ilev++)
			{
				auto dit_lev = Uconv[ilev].begin();
				auto dit_levp1 = Uconv[ilev + 1].begin();

				BoxData<double, 1> err = slice(Uconv[ilev][*dit_lev], varr);
				err -= Stencil<double>::AvgDown(2)(slice(Uconv[ilev + 1][*dit_levp1], varr));
				ErrMax[ilev] = err.absMax();
				std::string filename = "Comp_" + std::to_string(varr) + "_err_" + std::to_string(ilev);
				// NOTE: this assumes that the domain length is 1.0, which is assumed throughout this code. May cause errors if this changes.
				double dx = 1. / (err.box().size(0));
				if (inputs.saveConvTestData)
				{
					HDF5Handler h5;
					h5.writePatch({"err"}, 1, err, filename.c_str());
				}
				std::cout << "Lev: " << ilev << " , " << ErrMax[ilev] << std::endl;
			}
			double rate = log(abs(ErrMax[0] / ErrMax[1])) / log(2.0);
			std::cout << "order of accuracy = " << rate << std::endl;
		}
	}
	PR_TIMER_REPORT();
#ifdef PR_MPI
	MPI_Finalize();
#endif
}
