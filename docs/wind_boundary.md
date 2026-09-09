# Smooth super-fast wind boundary

This branch preserves the experimental boundary work. Its smooth analytic
verification does not establish fourth-order accuracy for production solar
wind. The supplied HDF5-driven SW study below did not show fourth order, and
an original-versus-new BC comparison has not yet established a production
advantage.

Enable the new radial boundary treatment for the file-driven solar wind
(`init_condition_type 0`) or analytic outflow (`1`) with:

```
-wind_boundary_order 4
-temporal_order 4
```

These are the only new settings needed with the usual `spatial_operator 1`
and `divergence_cleaning 1` configuration. No separate test executable or
wind-specific script is needed.

The input templates set `wind_boundary_order 4`; it is ignored for cases other
than 0, 1, 10, and 11. Omitting the option defaults to `0`, which retains the
existing treatment for cases 0 and 1. Cases 10 and 11 require `4`. Time order remains
a separate setting: use `temporal_order 4` for RK4. The new
implementation uses HOST data, radial coordinate 0, the existing exponential
radial mapping, at least ten angular and five radial cells per block, and
positive radial ghost-center radii. It checks
`u_r > c_fast,r` at the inner boundary. Being super-Alfvenic alone does not
satisfy this condition. RK3 and forward Euler are supported, but retain their
respective time orders. Use `limiter 0`, and inactive artificial
viscosity/floors for a smooth accuracy study. The existing `limiter 1` can
reduce maximum-norm convergence near smooth magnetic extrema.

## Implementation

`src/WindBoundary.H` contains the numerical kernels; `src/BoxOp_WindBoundary.H`
connects them to the solver. The boundary polynomials use dimensionless
conserved variables with radial momentum and two cubed-sphere tangential
momentum components, and the same decomposition of B as the interior solver.
Total energy remains a scalar; kinetic and magnetic energies use the physical
metric, since the tangential basis is not orthonormal.

The Powell MHD equations determine the boundary-normal derivative from the
prescribed time derivative, tangential derivatives, and gravity. A quartic
polynomial matches the boundary value, this normal derivative, and three
interior radial cell averages. Its integral supplies radial ghost averages;
angular averaging supplies the solver's cubed-sphere representation directly.
Numerical angular differentiation acts on local components. Analytic first
and second derivatives of the cubed-sphere direction supply the basis terms:
`d(T U)/deta = T dU/deta + (dT/deta) U`. `LocalFrame::rate` evaluates the
physical MHD RHS with these derivatives and rotates the resulting time
derivative back into the fixed local basis. This avoids numerically
differentiating the changing Cartesian direction of a purely radial flow.
Time derivatives, polynomial extrapolation, and RK stage evolution all use
the local components. No components are forced to zero during evolution.

The outer boundary already uses these local components and retains quartic
extrapolation of cell averages.
If extrapolation gives nonpositive density or internal energy, all ghost states
in that column are blended toward the last interior average. This leaves
admissible quartic extrapolations unchanged and reduces order where activated.
The solver reports these activations as `positivity-limited outflow column-stages`.
Auxiliary angular corners use a rotated nearby physical state if their anchor
is inadmissible.

If a reconstructed inner polynomial has nonpositive density or internal
energy, its variation is damped toward the boundary trace. This preserves the
trace, but locally reduces order when activated. The solver reports the number
of limited physical column-stages; the smooth accuracy test requires zero
activations there. Auxiliary angular ghost columns also need admissible
states because intermediate stencil calculations use rectangular arrays;
their fallbacks are counted separately. These do not impose extra physical
boundary conditions on the sphere.

The polynomial is advanced with the solver's actual RK stage formulas. Simply
evaluating the physical boundary function at each stage time would not, in
general, give consistent intermediate RK states. Gravity is included in each
RK stage and its mapped source is averaged to fourth order.

The new mode also corrects products of averaged quantities in the mapped
vector conversion, scalar fluxes, and Powell magnetic-divergence calculation.
The vector conversion adapter is local to this repository; it does not require
changes to the external Proto checkout.

For file data, spatial interpolation uses four samples in each angular
direction. Vectors are interpolated in Cartesian components, with periodic
longitude and pole reflection. Up to four frames are used for time
interpolation; fewer frames cannot provide fourth-order accuracy for arbitrary
time-varying data. Outside the supplied time range, the endpoint frame is
held while rotation continues. Five nearby evaluations supply the boundary
time derivative. Accuracy with respect to an underlying physical boundary
also requires adequate source-map resolution and frame cadence.

This treatment assumes smooth, compatible initial and boundary data. The
existing problem 1 initialization (`rho ~ r^-2`, `p ~ r^-2 gamma`, constant
radial speed) is not an exact steady solution with finite pressure. Its startup
is therefore not an exact-solution convergence test. Shocks, unresolved map
features, boundary time knots, and active limiters can reduce the observed
order.

## Accuracy test

Problems 10 and 11 share a smooth exact expansion and a density pattern
continuously entering through the inner boundary. Case 10 is hydrodynamic
(`B0=0`); case 11 has the previous nonzero, spatially uniform Cartesian field.
With `lambda = 1 + t/tau`, `tau = r0/(500 km/s)`, and initial wavelength
`L0 = (r1 - r0)/2`, it has

```
rho/rho0 = lambda^-3 [1 + 0.1 sin(2 pi (r/lambda - r0)/L0)]
v        = x/(tau + t)
p        = p0 lambda^(-3 gamma)
B        = B0 lambda^-2
```

There are two density cycles across the shell initially, independent of
`max_time`. The wavelength grows as `lambda L0`, leaving one cycle at `t = tau`.
Initialization, inner-boundary data, and exact-error diagnostics all call the
same `WindBoundary::exact` function.

It solves ideal MHD without gravity: pressure and magnetic field are spatially
uniform, and the velocity has zero material acceleration. Case 11 evolves
the full magnetic field and subtracts the numerical momentum-force
residual of a uniform magnetic equilibrium, as described below. Case 10
has no magnetic field or equilibrium subtraction. Errors
compare the evolved mapped cell averages against an
independent tensor Gauss quadrature of the exact solution, including the cells
next to the inner boundary.

Use the existing convergence scripts and `exec/inputs_convergence`. Both the
local and Anvil scripts accept problems 0, 1, 10, and 11 alongside the pulse cases
2–5. They preserve `wind_boundary_order` and `temporal_order` from the input
file; the existing `HELIOCUBED_CONVERGENCE_TEMPORAL_ORDER` override is optional.

For a smooth wind accuracy study, set these values in `exec/inputs_convergence`:

```
-wind_boundary_order 4
-temporal_order 4
-limiter 0
-Sun_gravity 0
-CFL 0.4
-max_time 1496.4
```

Then select problem 10 (zero field) or 11 (uniform field) in the ordinary
convergence driver:

```sh
HELIOCUBED_CONVERGENCE_PROBLEM_TYPE=10 \
HELIOCUBED_CONVERGENCE_TEST_TYPE=2 \
HELIOCUBED_CONVERGENCE_MAX_ITER=1000 \
bash scripts/run_convergence_test.sh
```

Set the mesh, box sizes, and MPI counts in that script's existing configuration
as usual. Mode 2 refines space and time together; mode 1 holds the timestep
fixed. Use enough iterations to reach the same `max_time` at all resolutions.
The executable uses the shell geometry from its normal build; coarse radial
grids must still have positive ghost-center radii.

On Anvil, select `TEST_CASES=(10)` (or
`HELIOCUBED_CONVERGENCE_TEST_CASES=10`) in
`scripts/run_convergence_test_anvil.sh` and set its existing mesh, rank, and
per-case time settings. The default problem-10 time limit is 1496.4 seconds.

The common driver writes the usual three-resolution comparison to
`Convergence_results/exact_wind/convergence_summary.txt`, together with each
level's exact L1, Linf, and inner-boundary Linf errors. Anvil writes the same
files beneath its job-specific results directory. Check that physical boundary
positivity limiting stays inactive when measuring smooth convergence.

Select problem 1 for analytic-outflow self-convergence, or problem 0 for
file-driven wind with an absolute `-BC_file` path in the input. Use the desired
gravity setting for those cases. Their startup compatibility and map
smoothness affect the measured order, as discussed above.

## Running on the full solar-wind shell

`exec/inputs_smooth_outflow` selects case 10 (zero field);
`exec/inputs_smooth_outflow_const_B` selects case 11 (uniform Cartesian field).
Both use RK4, `CFL 0.4`, no artificial
viscosity or reconstruction limiter, a 30-by-90 mesh, and a final time of
29928 seconds. On the 21.5–231.5 solar-radius shell, the density wavelength
starts at 105 solar radii (two cycles) and ends at 210 solar radii (one cycle).
Both inputs use `scale_data_with_r 0` to show the physical density in g/cm^3.
The solar-wind output option `scale_data_with_r 1` instead writes
`(rho/m_p) (r/AU)^2`; that radial factor masks this test's 10% density wave.
Changing the option affects newly written output, not existing HDF5 files.
Refine the radial grid for an accuracy study, and check both inflow and outflow
limiting counts.

For a two-to-one-cycle run through the local convergence driver, select the
corresponding file with `HELIOCUBED_CONVERGENCE_INPUT`, set
`HELIOCUBED_CONVERGENCE_PROBLEM_TYPE` to 10 or 11, and allow enough iterations
with `HELIOCUBED_CONVERGENCE_MAX_ITER` (for example 1000). The driver’s default
three-step limit is only a short check. Anvil now defaults both of these
cases to 29928 s. If the compiled inner radius changes, use
`max_time = r0/(500 km/s)` to retain the two-to-one-cycle interval.

The case split and short regression checks are recorded in the
[case-10/11 results](../Convergence_results/problem10_11_split_20260909/README.md).

The original failure on a 15-by-45 mesh with the previous, shorter initial
wavelength (`L0 = r0`) came from negative density in quartic outer
ghost extrapolation; changing Euler to RK4 alone did not fix it. The outflow
positivity treatment above prevents that extrapolation failure. The timestep
check also stops a run if its timestep becomes nonfinite, nonpositive, or too
small to advance time.

After pulling these changes, rebuild once with
`make -C exec -f GNUmakefile_TS -B -j4` (use the corresponding makefile on your
cluster). Header dependency targets were corrected so subsequent header edits
trigger normal incremental rebuilds.

## Recorded validation

The implementation was validated on 2026-09-09 with six local MPI ranks, RK4,
`CFL 0.4`, and `limiter 0`, using 16, 32, and 64 cells per direction on all
six faces. That study used a short shell (`r0 = 1496400000000 cm`, `r1 = 2 r0`)
and final time 1496.4 seconds. It used the previous density wavelength
`L0 = r0`, not the current two-cycle initialization. These recorded results
are specific to that geometry and profile. The measured orders from 32 to
64 cells per direction are below.
The inner norm covers the first three radial layers.

| Conserved variable | Volume-weighted L1 | Global Linf | Inner Linf |
| --- | ---: | ---: | ---: |
| Density | 3.769 | 4.291 | 4.715 |
| Cartesian momentum x | 3.610 | 4.316 | 3.791 |
| Cartesian momentum y | 3.610 | 4.316 | 3.792 |
| Cartesian momentum z | 3.610 | 4.316 | 3.792 |
| Total energy | 3.654 | 4.318 | 4.011 |
| Cartesian Bx | 4.411 | 4.023 | 4.197 |
| Cartesian By | 4.432 | 4.066 | 4.251 |
| Cartesian Bz | 4.445 | 4.140 | 4.241 |

No physical inflow columns required positivity limiting in this accuracy
study. Auxiliary angular ghost fallbacks occurred on the coarse grids and
were inactive at 64 cells. The same 16-cell physical mesh split into smaller
patches passed the decomposition comparison.

Additional checks passed: manufactured MHD equations and gravitational work;
super-fast/reversed-flow checks; ghost-average convergence; RK3/RK4 local
errors; boundary-image interpolation through longitude seams and poles;
mapped vector conversion including the radial end cells; and the five
existing time-limit/output regression tests against the rebuilt production
executable.

Problem 1 completed three steps with gravity under both RK3 and RK4. The
supplied `psi_eclipse24_swig_wsa2_smth5_bc_r21_5_corotating_Bp_two_frames.h5`
completed five RK4 steps with gravity on a 32-cell-per-direction short shell,
reaching approximately 244.2 seconds. These runs passed the finite-state
check and activated positivity limiting. They are startup/robustness checks,
not long-duration production validation or fourth-order studies of that map.

The original validation artifacts remain locally in `Test_results/wind_boundary/`
(`convergence.json` and run inputs/logs). The standalone test sources and runner
were subsequently removed; ongoing wind convergence runs use the shared
convergence scripts above.

## Validation after the outflow fix

The updated `inputs_smooth_outflow` completed with six MPI ranks, including its
normal output and probe settings: 11 RK4 steps reached 1496.4 seconds without
timestep collapse. Short runs on 30-by-90 and 30-by-180 meshes also completed;
outflow limiting remained active, so these are not fourth-order convergence
measurements.

The resolved 32-cell short-shell benchmark reproduced all 24 previously
recorded exact-error norms unchanged, with zero inflow and outflow limiting.
The 16-cell short-shell benchmark now activates outflow limiting and its errors
change; the table above records the earlier 32-to-64 study. The five-step
file-driven SW check and all five existing time-limit regressions also passed.

## Historical validation

The reports below predate the case split: their magnetized "case 10" is now
case 11. Those archived inputs need `-init_condition_type 11` to reproduce
the magnetic test with the current executable. The same applies when
continuing an older magnetized case-10 checkpoint.

## Validation of radial-component boundary update

The 2026-09-09 radial-component refactor was checked with analytic geometry
derivatives, a general Cartesian linear-state MHD RHS including gravity, and
a radial hydrodynamic ILW slope. The slope agrees with the analytic value to
roundoff. In a full-shell 30-by-90 zero-field case 10 control, the maximum
transverse speed after one step decreased from about 0.054 cm/s to
3.8e-7 cm/s. Short RK4 runs for cases 0, 1, and 10 also completed.

The existing three-level convergence mode was run on the compact shell
`r1 = 2 r0`, with the current two-cycle profile, using 16-by-64, 32-by-128,
and 64-by-256 meshes per block. All levels reached 31.7170677 seconds with
the same nominal timestep. No physical boundary positivity limiting was
active. The finer-pair L1 orders range from 3.92 to 4.15 across all eight
conserved variables; global Linf orders range from 3.84 to 3.89. This is a
short-time spatial check, not a long-time solar-wind accuracy guarantee.

The magnetized case 10 still has small transverse errors: its full-shell
one-step maximum was 0.155 cm/s, versus 0.140 cm/s before this change.
The initial pressure reconstruction error is also unchanged. Radial
hydrodynamic symmetry and exact preservation of a constant Cartesian B are
different discrete properties; the latter still involves approximate
magnetic/metric cancellations.

Inputs, run logs, error tables, and limitations are recorded in
`Convergence_results/problem10_radial_bc_20260909/README.md`.

## Full-shell magnetic seam audit

The subsequent full-shell audit reproduced the reported case-10 seams:
maximum `abs(Vp)` was 3.097 cm/s at 5435.24 s on the 30-by-90 mesh.
A separate projection of the solver state confirmed that this is primarily
an evolved-state error, rather than a plotting conversion error. Supplying
analytic local cell averages in a diagnostic step still generated transverse
velocity, implicating the mapped magnetic operator as well as its interface
treatment. It is not explained by a radial/tangential BC variable choice alone.

At 119.244232 s, refining 30-by-90 to 60-by-180 with nearly matched timesteps
gave observed orders of 3.16 for maximum `abs(Vp)` and 3.30 for maximum
transverse speed. This two-mesh comparison does not establish asymptotic
orders, but it does not demonstrate uniform fourth-order accuracy of these
local diagnostics. It does not invalidate the earlier compact conservative
L1 study, which used different norms and a different shell.

Wider BC derivatives, Cartesian magnetic exchange, higher-order inter-block
interpolation, and equilibrium subtraction were tried in isolated builds.
None was established as a general cure; no trial solver change was retained
at that stage of the investigation.
In particular, the existing constant-B cases 3, 5, 7, and 9 already subtract
a numerical magnetic-equilibrium residual, whereas that case-10 audit used
the full-field operator without that correction. Their clean symmetry is not
evidence that the ordinary operator exactly preserves a uniform Cartesian B.

See [the seam audit](../Convergence_results/problem10_seam_audit_20260909/README.md)
for measurements, unsuccessful correction attempts, and remaining work.

## Case-11 equilibrium force subtraction

Following that audit, case 11 now automatically uses equilibrium force
subtraction. No extra input option is required. Cases 0, 1, and 10 are unaffected;
the existing constant-field pulse cases retain their previous subtraction.

`WindBoundary::magneticEquilibriumAverage` constructs angular cell averages
of a stationary, spatially uniform Cartesian magnetic reference, with zero
velocity and uniform positive pressure. The same mapped operator evaluates
its numerical residual, which is cached per patch/operator instance. Since
this reference has zero continuum force, subtracting its discrete momentum
residual removes part of the geometry error. The corresponding momentum
fluxes are also subtracted before interface flux registration, so refluxing
uses the corrected fluxes.

The case-11 field is `B(t)=B0/lambda^2`. The cached reference force and
momentum fluxes are therefore multiplied by `lambda^-4`, using the actual
RK stage time (also on a restart). The stationary reference pressure is
scaled with its magnetic energy for this force calculation; it is not a
replacement for the expanding solution's physical pressure.

Case 11 continues to evolve the full state. No velocity component is reset,
and its density, energy, and induction equations receive no direct reference
subtraction. This also avoids subtracting any numerical induction/energy
diffusion that a reconstruction limiter might produce for the reference.
Boundary values, initialization, and independent exact-error quadrature
are unchanged.

This is a case-specific improvement, not a general SW correction or exact
preservation of radial symmetry. Residual seams and pressure reconstruction
errors remain. Convergence results obtained with this correction should be
identified as using equilibrium subtraction.

Validation is recorded in
[the equilibrium-subtraction results](../Convergence_results/problem10_equilibrium_20260909/README.md).

## HDF5-driven SW study

The two-frame file in `scripts/` was tested on the full 21.5–231.5 solar-radius
shell with RK4 and boundary order 4. The three meshes were 15×15×45, 30×30×90,
and 60×60×180 per block, all advanced to 300 seconds with a common timestep.
Global volume-weighted L1 self-convergence orders, including both boundaries:

| Conserved variable | SW limiter/viscosity settings | Limiter/viscosity disabled |
|---|---:|---:|
| Density | 1.57 | 1.55 |
| Momentum components | 1.50–1.62 | 1.47–1.58 |
| Total energy | 0.41 | 0.38 |
| Magnetic components | 1.64–1.79 | 1.65–1.79 |

Low order was already present at initialization. The input contains sharp
magnetic polarity reversals and other poorly resolved structure. Boundary
positivity limiting remained active at all resolutions in both configurations.
The coarse primitive output also had negative pressures, although the
conservative mean states remained admissible; medium and fine final primitive
outputs had positive pressure. Smaller-timestep controls did not materially
change the coarse-pair L1 refinement differences.

This was a short startup study, not a test of a relaxed wind or an exact-solution
error measurement. Both configurations used the new BC, so the study does not
compare its accuracy against the original BC. Full local inputs, logs, saved
states, and analysis are in `Convergence_results/sw_h5_accuracy_20260909/`;
generated results are excluded from Git.

## References

- Tan and Shu, *Inverse Lax-Wendroff procedure for numerical boundary conditions
  of conservation laws*, JCP 229 (2010), 8144–8166.
- Zhao, Huang and Ruuth, *Boundary treatment of high order Runge-Kutta methods
  for hyperbolic conservation laws*, JCP 421 (2020), 109697;
  [arXiv:2001.09854](https://arxiv.org/abs/2001.09854).
