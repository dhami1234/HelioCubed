# HelioCubed

HelioCubed solves the ideal magnetohydrodynamic (MHD) equations on a six-panel
cubed-sphere shell using the [Proto](https://github.com/applied-numerical-algorithms-group-lbnl/proto)
framework. It supports solar-wind and coronal mass ejection (CME) simulations,
idealized test problems, and convergence studies, with MPI parallelism and
adaptive mesh refinement (AMR).

The `CubedSphere` branch includes both uniform-grid and AMR solvers:

- Second-order piecewise linear reconstruction with minmod or MC limiting.
- Fourth-order finite-volume reconstruction with the standard mapped scheme
  ([McCorquodale et al., 2015](https://doi.org/10.1016/j.jcp.2015.01.006)) or
  pointwise Berta WENOZ reconstruction
  ([Berta et al., 2024](https://arxiv.org/abs/2310.11831)).
- Rusanov and Roe eight-wave Riemann solvers.
- Multilevel AMR with conservative prolongation, restriction and refluxing;
  separate radial and angular refinement ratios; optional time subcycling.
- File-driven or analytic solar wind, solar gravity, scheduled FRiED CME insertion,
  trajectory probes, slices, and mapped HDF5 output.

For high-order refinement background, see
[McCorquodale and Colella (2011)](https://doi.org/10.2140/camcos.2011.6.1),
*A high-order finite-volume method for conservation laws on locally refined grids*,
and [McCorquodale, Ullrich, Johansen, and Colella (2015)](https://doi.org/10.2140/camcos.2015.10.121)
on adaptive multiblock methods for the shallow-water equations on the sphere.

## Build

You need a C++ compiler compatible with the selected makefile (C++20 for the
default and local builds), GNU Make, Proto, HDF5, and BLAS/LAPACK. MPI builds also
need an MPI compiler wrapper and parallel HDF5 built with the same MPI stack.

Use Proto's `integrate/master-cubedSphereDev` branch, which includes the
multiblock AMR, mapped HDF5, and checkpoint support used here. Proto is a separate
checkout; it is not downloaded by the HelioCubed makefiles.

```sh
git clone --branch CubedSphere https://github.com/dhami1234/HelioCubed.git
git clone --branch integrate/master-cubedSphereDev \
  https://github.com/applied-numerical-algorithms-group-lbnl/proto.git
cd HelioCubed/exec
```

Choose a makefile and adjust its compiler and dependency paths for your system:

| Makefile in `exec/` | Configuration |
| --- | --- |
| `GNUmakefile` | Default Clang, serial build |
| `GNUmakefile.clang` | Alternate Clang, serial build |
| `GNUmakefile.gcc` | GCC, serial build |
| `GNUmakefile_TS` | Local macOS MPI build with Homebrew HDF5/OpenBLAS |
| `GNUmakefile_Anvil` | Anvil MPI build |
| `GNUmakefile_Pleiades` | Pleiades MPI build |

For example, build the default serial executable with explicit library paths:

```sh
make -B -j8 PROTO=/absolute/path/to/proto \
  HDF5_INC=/path/to/hdf5/include HDF5_LIB=/path/to/hdf5/lib \
  LAPACK_LIB="-L/path/to/openblas/lib -lopenblas"
```

For a configured MPI makefile, use, for example:

```sh
make -f GNUmakefile_TS -B -j8
```

The executable is `exec/cubedSphereTest.exe`. All supplied makefiles enable
`PR_AMR`; select uniform-grid or adaptive operation at runtime with
`amr_max_level`. Rebuild with `-B` after changing compiler flags or makefiles so
existing objects are recompiled. The makefiles are the supported standalone
build path; the CMake fragments are for integration with a parent build.

## First run

Run from `exec/`. Input files use `-option value`, with comments on separate
lines. Paths in inputs are relative to the directory where the executable runs.
The program accepts either one input filename or command-line options; it does
not combine an input filename with extra command-line overrides.

Create an `inputs_run` file with this short analytic solar-wind example, which
does not need boundary data or a checkpoint:

```text
-init_condition_type       1
-strength                  100
-Sun_gravity               1
-gamma                     1.5
-domainSize                24
-thickness                 24
-inner_radius              1.4964e12
-outer_radius              1.61124e13
-radial_stretch            0
-boxSize_nonrad            12
-boxSize_rad               12
-spatial_order             2
-limiter                   3
-Riemann_solver_type       2
-temporal_order            2
-linear_viscosity          0
-nonlinear_viscosity       0
-boundary_reflux           1
-CFL                       0.5
-max_iter                  5
-max_time                  3600
-amr_max_level             0
-P_floor_cadence           1
-write_cadence             1000
-write_time_cadence        1000000000
-checkpoint_cadence        1000
-data_file_prefix          quickstart
-checkpoint_file_prefix    quickstart_checkpoint
```

```sh
./cubedSphereTest.exe inputs_run                 # serial build
mpirun -np 8 ./cubedSphereTest.exe inputs_run     # MPI build
```

Volume output is written at initialization and at the final step. `max_iter`
and `max_time` are absolute stopping limits, including elapsed steps/time when
restarting. Time inputs are in seconds and shell radii are in centimeters.
Every input file must explicitly set `inner_radius`, `outer_radius`, and
`radial_stretch`; a missing or commented setting stops the run with an error.
Use `radial_stretch=0` for uniform radial spacing or a nonzero value for
exponential stretching. These settings do not come from the Makefile.

The established input files provide more complete configurations:

| Input file | Purpose |
| --- | --- |
| [inputs_SW](exec/inputs_SW) | File-driven solar wind, CME, probes, and AMR |
| [inputs_Sine_Wave](exec/inputs_Sine_Wave) | Analytic solar wind |
| [inputs_Test_Problems](exec/inputs_Test_Problems) | Pulse, shock, and blast problems |
| [inputs_convergence](exec/inputs_convergence) | Spatial and space/time convergence studies |
| [inputs_local](exec/inputs_local) | Local run configuration |

`inputs_SW` starts from initial conditions by default and schedules a CME;
its example `restart_file` setting is commented out. Before running, check the
grid, physical epoch, boundary file, CME settings, and output prefixes.
Uncomment `restart_file` only when the matching checkpoint set is available.
With `init_condition_type=0`,
set `BC_file` to a compatible boundary HDF5 file. The repository includes one
example in `scripts/` and a [SWiG conversion script](scripts/SWiG_to_HelioCubed_h5_conversion.py).
`BC_start_time` and `CME_Enter_Time` accept decimal years or
`DD-MM-YYYYTHH:MM:SS` in UTC.

## Adaptive mesh refinement

`amr_max_level=0` uses the uniform-grid driver. A positive value enables AMR:
`1` allows levels L0 and L1, `2` allows L0 through L2, and so on.

| Spatial method | AMR time integrators | Required settings |
| --- | --- | --- |
| O2 | `temporal_order=2` or `3` | `boundary_reflux=1`; artificial viscosity off |
| O4 standard or Berta | `temporal_order=3` or `4` | `boundary_reflux=1`, `state_conversion=1` |

AMR base dimensions must be at least 12 cells and divisible by their patch
sizes, which must be at least 8. O4 additionally requires at least 20 angular
cells per panel. `domainSize` is the angular cell count **on each panel**;
`thickness` is the radial cell count. The base grid has
`6 * domainSize * domainSize * thickness` cells.

To make the first-run example adaptive, change its `amr_max_level` to `1` and
add the following settings:

```text
-amr_ref_ratio_rad          2
-amr_ref_ratio_nonrad       1
-amr_subcycling            1
-amr_regrid_interval       4
-amr_regrid_time_cadence    3600
-amr_buffer                2
-amr_verbosity             1
-amr_refine_density_thresh -1
-amr_refine_pressure_thresh -1
-amr_refine_B_thresh       -1
-amr_refine_V_thresh       0.2
```

Each ratio can be one integer for every interface or a comma-separated list
with one entry per interface. For example, `amr_max_level=2`,
`amr_ref_ratio_rad=2,4`, and `amr_ref_ratio_nonrad=1,2` refine radially by 2
then 4, and angularly by 1 then 2. Ratios such as 2, 4, and 8 are supported;
1 leaves a direction unchanged. Every interface must refine at least one
direction, and ratios must allow aligned fine patches to tile the domain.

Fine patches can have different cell counts from base patches. Set
`amr_boxSize_rad` and `amr_boxSize_nonrad` to one value for every fine level
or a comma-separated list. Zero retains the existing automatic patch sizing.
Positive sizes must be at least 8 cells, divisible by that level's refinement
ratio, and tile both the refined domain and a refined parent patch; a refined
patch must contain at least two parent cells in each refined direction.

For a base grid of 45 angular by 120 radial cells, base patches of 15 angular
by 30 radial cells, and refinement ratios 2 angular / 4 radial, use:

```text
-amr_boxSize_rad            20
-amr_boxSize_nonrad         0
```

This divides a 60×30×30 fine patch into three patches of 20×30×30. Setting
`amr_boxSize_nonrad=10` instead creates 27 patches of 20×10×10 and spreads very
sparse fine work across more ranks. Smaller patches increase ghost-cell and
message overhead, so the best sizes depend on the number of ranks and fine
patches. Automatic tagging may cover a different region when the patch
dimensions change; fixed regions that tile both layouts allow a comparison
with identical refined cells. Patch dimensions are printed in the AMR grid
summary.

With `amr_subcycling=0`, all levels share the finest required timestep. With
`1`, finer levels take extra steps only when their CFL limits require them.
RK3/RK4 subcycling uses extrapolated SSPRK2 advances, which adds work; its
performance depends on the refined volume and the level timestep limits.

Refinement uses the MS-FLUKSS relative-gradient magnitude of density, pressure,
magnetic-field magnitude, or speed. For an interior cell, each directional
contribution is `(q_plus - q_minus) / (q_plus + q_minus)`; the criterion uses
the Euclidean norm of those contributions. Physical radial boundaries use
one-sided differences. Thus, a threshold of `0.2` is a relative-gradient test,
not a direct 20% difference from each adjacent cell.

Any enabled criterion can tag a cell. Negative thresholds disable a criterion;
zero tags any nonzero gradient. The corresponding
`amr_refine_<variable>_maxlevel` options (`density`, `pressure`, `B`, `V`) limit
the finest child level: `-1` uses `amr_max_level`, and `0` disables that
criterion. Tags grow by `amr_buffer` parent cells and allocate whole fine
patches. Large patches can therefore refine much more volume than the tagged
cells alone suggest.

Optional physical bounds restrict where cells can trigger refinement:

```text
-amr_refine_r_range       1.4964e12,7.482e12
-amr_refine_colat_range   60,120
-amr_refine_lon_range     350,20
```

Radius is in cm. Colatitude is in degrees from the north pole: 0 north,
90 equator, 180 south. Longitude is in degrees in the simulation's Cartesian
frame: 0 along +X, 90 along +Y. These are fixed directions in that frame;
boundary-frame rotation does not rotate the bounds. Longitudes are normalized
modulo 360, so negative values work and `350,20` crosses 0. A span of at least
360 degrees, such as `0,360` or `-180,180`, allows every longitude; equal
endpoints select a single meridian. Longitude does not restrict either pole.

Omit any range to leave that coordinate unrestricted. All supplied ranges
must contain a parent cell's center for it to tag. The bounds apply to every
AMR level and to both variable-dependent criteria and `amr_refine_regions`;
they restrict refinement without forcing it. Buffering, whole-patch allocation
and nesting support can extend the actual refined mesh beyond the bounds.
Changing bounds on restart takes effect at the next regrid; the saved mesh
is restored first. Commented examples are included in the established inputs.

Regridding occurs when either cadence is reached: `amr_regrid_interval` counts
base steps, and `amr_regrid_time_cadence` measures elapsed simulation seconds.
Zero disables that cadence. Time-based regrids occur at the first synchronized
coarse step crossing a cadence boundary. Verbosity 1 prints level timesteps,
regrid duration, and coverage as percentages of base-grid volume, including
regrid attempts that leave the hierarchy unchanged; 2 adds detailed advance
and reflux events, while 0 retains ordinary iteration progress.

`P_floor_cadence` also controls AMR density repair. Density uses the existing
`RHO00` floor and pressure has a zero lower bound. Both are checked before
regridding so newly refined parents are admissible; repairs are reported when
they occur.

## Output, VisIt, and restart

AMR volume output is a native multilevel HDF5 hierarchy:

```text
data_amr_0.hdf5
data_amr_0.map.hdf5
```

Keep each data file beside its matching map file. Open the **data** file in
VisIt using its Chombo reader. In the reader options, enable checking for a
mapping file and disable enabling only the root level. For the physical shell,
apply a Displace operator with `_mapping_displacement` and factor 1 to each
plot; the reader creates that expression from the map file. These settings are
supported by the [VisIt Chombo reader](https://github.com/visit-dav/visit/blob/develop/src/databases/Chombo/avtChomboFileFormat.C).

Plot `amr_level` to see which levels exist; `active=1` identifies cells not
covered by finer data. To inspect the refined mesh, enable all available levels
in the Subset/SIL selection and raise any MultiresControl resolution setting to
the finest level before adding a Mesh plot. A file may contain only L0 if no
criterion has tagged cells.

`output_coordinates` selects the vector components written in **both volume
files and slices**, with the same behavior for uniform grids and AMR:

| Setting | Velocity fields | Magnetic fields |
| --- | --- | --- |
| `0` (default) | `Vr`, `Vt`, `Vp` | `Br`, `Bt`, `Bp` |
| `1` | `Vx`, `Vy`, `Vz` | `Bx`, `By`, `Bz` |

Spherical components use the orthonormal radial, colatitude, and azimuthal
directions. Both choices also write `density` and `P`; AMR volumes add `active`
and `amr_level`. This setting controls output independently of the solver's
`state_conversion` setting. The mesh geometry remains the physical shell.

`scale_data_with_r=0` writes raw CGS values in either coordinate basis. With
`1`, write density multiplied by `(r/AU)^2/c_MP`, pressure by
`(r/AU)^(2*gamma)*1e12`, and velocity in km/s. Magnetic fields are converted to
microgauss with radial components multiplied by `(r/AU)^2` and transverse
components by `r/AU`. For Cartesian output this magnetic scaling is applied
in the spherical basis before rotating back to Cartesian components.
`output_divB=1` always adds physical `divB` in G/cm.

`slices` selects any combination of `X Y Z R`; AMR slices omit covered coarse
cells. Checkpoints retain conservative state. Probes sample the finest available
level at each trajectory position. The bundled solar-wind slice plotting
defaults use spherical field names.

Probes write orthonormal RTN velocity components (`VR`, `VT`, `VN`) and magnetic
components (`BR`, `BT`, `BN`) for both uniform-grid and AMR runs, independently
of the volume/slice `output_coordinates` setting. For plotting in spherical
components, use `Vr=VR`, `Vtheta=-VN`, and `Vphi=VT`, with the same mapping for B.
Theta is colatitude, increasing southward.
Each headerless probe row contains decimal year, radius (AU), latitude (degrees),
simulation longitude (degrees), number density (cm^-3), three velocity components
(km/s), pressure (picodyne/cm^2), and three magnetic components (microgauss).
Probes do not apply the volume/slice radial scaling.

Slices from both uniform-grid and AMR runs are self-contained **Tecplot ASCII
`.dat` files**, for example `data.z.100.dat` or `data.r.100.dat`. They include
physical coordinates, mesh connectivity, and cell-centered primitive fields.
`X`, `Y`, and `Z` select the planes x=0, y=0, and z=0; `R` selects the spherical
shell at radius `1.5e13` cm. Slice cadence is controlled by `slice_cadence` and
`slice_time_cadence`.

- **Tecplot:** use **File > Load Data**, choose the Tecplot data loader, and
  open the slice `.dat` file. Enable contours for a field such as `density` or
  `P`, and mesh lines to inspect the refinement. See the
  [Tecplot loading guide](https://tecplot.azureedge.net/products/360/current/help/users_manual/data-loading.html).
- **ParaView:** use **File > Open**, select the slice file or file series,
  choose **Tecplot Reader** if prompted, and click **Apply**. Color by a cell
  field and use **Surface With Edges** to inspect the mesh. ParaView's
  [Tecplot Reader](https://www.paraview.org/paraview-docs/v5.13.0/python/paraview.simple.TecplotReader.html)
  supports multiple zones, block-packed cell data, and file series.

An AMR restart checkpoint is a **set** of files:

```text
checkpoint_amr_100.amr
checkpoint_amr_100_L0.hdf5
checkpoint_amr_100_L1.hdf5
```

Keep the manifest and every matching `_L*.hdf5` file together, and set:

```text
-restart_file checkpoint_amr_100.amr
```

Keep the original base grid, patch sizes, shell geometry, gamma,
`init_condition_type`, configured maximum AMR level, and refinement ratios.
New checkpoints also save the fine patch sizes. Zero-valued fine-size inputs
restore those saved dimensions; explicit positive sizes must match them.
For checkpoints written before fine-size support, keep both fine-size inputs
at zero to retain the original automatic layout.
The restart restores the saved hierarchy, iteration, simulation time, and CME
insertion state. Set stopping limits beyond the saved step/time to continue.
Plot files are for visualization; restart uses the checkpoint manifest and
its conservative-state files. `checkpoint_cadence` controls checkpoint output,
and `max_checkpoint_files` controls retention of sets written during the run.

## Convergence and analysis

The same convergence driver supports uniform grids and AMR:

- `convTestType=1`: refine space while holding the timestep fixed.
- `convTestType=2`: refine space and time together.
- `convTestType=3`: compare initialized states without advancing.
- `convTestType=4`: initialize one grid without advancing.

For modes 1/2, `convergence_level=-1` runs all three study resolutions;
`0`, `1`, or `2` runs one study, and `3` compares saved studies.
`convergence_dt` sets a positive fixed timestep or uses `-1` for automatic
selection. For split mode-2 runs, halve the prescribed timestep and increase
the step limit for each successive study to reach the same physical time.
AMR studies compare conservative composite projections onto each study's
base mesh. Study indices are separate from AMR levels inside each hierarchy.

[run_convergence_test.sh](scripts/run_convergence_test.sh) automates split
studies, and [run_test_problems.sh](scripts/run_test_problems.sh) runs idealized
cases and plots slices. Edit their configuration sections or documented
environment overrides before launching; the test-problem script replaces
existing result folders for the selected cases. With AMR, choose a supported
time integrator in the script settings and the required O4 conversion setting
in the input template. Their `_anvil.sh` variants handle batch runs.

Python plotting dependencies are listed in
[requirements-plotting.txt](scripts/requirements-plotting.txt). Timing analysis
scripts in `scripts/` process Proto `.time.table` reports; the AMR report is
`CubedSphere_AMR.time.table`.
