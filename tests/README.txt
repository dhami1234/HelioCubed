Boundary average regression

The file boundary path must provide native conserved CELL AVERAGES to
PreStagePatch, while initialization needs classical primitive POINT SAMPLES
at cell centers. boundary_averages.cpp tests both contracts without an SW file.

The test checks all six blocks, radial/angular corner ghosts, radial patch
source strips, longitude wraparound, and both node- and center-origin source
longitude grids. It compares uniform radial conserved averages to analytic
radial integrals, requires zero transverse momentum and magnetic components
for purely radial input, and checks centered interpolation for all eight
source variables against independently mapped cell centers. Nonfinite values
and scaled errors above 5e-12 fail even when NDEBUG is set.

Use your normal platform makefile and compiler/library configuration. From
the HelioCubed root, override its source and executable as follows:

  make -B -C exec -f GNUmakefile_Pleiades \
    TARGET=../tests/boundary_averages.exe \
    SRC=../tests/boundary_averages.cpp \
    VPATH='../tests:../src:/nobackup/tsingh4/proto/examples/_common'
  mpirun -np 12 ./tests/boundary_averages.exe tests/inputs_boundary_averages

Use the MPI launcher required by your machine's scheduler. Adjust VPATH to
match your Proto path if needed. The example builds a separate executable
and leaves cubedSphereTest.exe untouched. A serial build can be run directly.
The input fixture requires gamma=1.5 and radial_stretch=0 for its analytic
radial-integral check; boundary production code also supports other values.

The broader seam investigation separately compares initial numerical MHD
rates against an analytic physical RHS and evolves the real rotating SW
map. The regression above checks boundary-state correctness; by itself it
does not establish that every evolved seam artifact has disappeared.

Seam pressure regression

seam_pressure.cpp initializes a sharp radial velocity/Br front with known,
constant thermal pressure. It calls the production prepareState routine and
compares its pressure with that analytic value before any time step. It also
checks that density, momenta and magnetic fields everywhere, owned energy,
and physical radial ghosts (including angular corners) remain unchanged.
It requires energy corrections at both the first and last physical radial
layers, so a correction accidentally limited to the inner boundary fails.

Build as above with TARGET=../tests/seam_pressure.exe and
SRC=../tests/seam_pressure.cpp, then run:

  mpirun -np 12 ./tests/seam_pressure.exe tests/inputs_seam_pressure

Repeat with boxSize_nonrad=24, or a different rank count, to exercise patch
cuts through seams. The standard fixture uses 48 angular cells per block.
Failures return a nonzero exit status even when NDEBUG is set.

The seam_pressure_correction input defaults to 1. Setting it to 0 restores
the previous state preparation for comparison. The enabled path converts JU,
fills radial BCs, transfers recovered pressure with the existing interpolation
stencil, and adjusts only angular ghost energy within the radial domain. A
temporary extra halo makes the entire ordinary ghost region available for
pressure recovery; three corrections with fixed targets handle transitions
out of the existing positivity fallback. Nonfinite or nonpositive transferred
pressure is not used to change energy. The interpolation plan is cached while
temporary state buffers are released after preparation. Stepping, output
pressure recovery and CME insertion use the same preparation routine.

This changes interblock energy transfer. The pressure improvement does not
imply smaller errors in every conserved variable: in the rotating SW control,
ghost energy RMS error increased from 0.218% to 0.397% as pressure RMS error
fell from 2.07% to 0.111%. Velocity and magnetic interpolation errors remain.
The regression does not establish long-term stability or removal of evolved
Btheta seam structures, nor guarantee positivity for arbitrary shocks.
