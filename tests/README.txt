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
