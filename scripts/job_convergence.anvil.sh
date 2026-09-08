#!/bin/bash
#SBATCH --account=mca07s033
#SBATCH --partition=wholenode
#SBATCH --nodes=9
#SBATCH --ntasks=1152
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --job-name=hc-conv-B
#SBATCH --output=hc-conv-B.o%j
#SBATCH --error=hc-conv-B.e%j

set -Eeuo pipefail

trap 'rc=$?; echo "ERROR: line ${LINENO}: ${BASH_COMMAND} (exit ${rc})" >&2; exit "$rc"' ERR

# ---------------------------------------------------------------------------
# Software environment: compiler before MPI.
# ---------------------------------------------------------------------------

module purge
module load xalt/2.10.45
module load intel/19.0.5.281
module load numactl/2.0.14
module load libszip/2.1.1
module load zlib/1.2.11
module load mvapich2/2.3.6
module load hdf5/1.10.7

module list
command -v mpirun

export MV2_HOMOGENEOUS_CLUSTER=1
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# ---------------------------------------------------------------------------
# Shared-library paths.
# ---------------------------------------------------------------------------

export LD_LIBRARY_PATH="/apps/spack/anvil/apps/hdf5/1.10.7-intel-19.0.5-gddu73a/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export LD_LIBRARY_PATH="/home/x-dhami/gsl/lib:${LD_LIBRARY_PATH}"

# OpenBLAS location specified in exec/GNUmakefile_Anvil.
OPENBLAS_LIB="/apps/spack/anvil/apps/openblas/0.3.17-gcc-11.2.0-2qrsari/lib"

if [[ ! -r "${OPENBLAS_LIB}/libopenblas.so.0" ]]; then
    echo "ERROR: Cannot find ${OPENBLAS_LIB}/libopenblas.so.0" >&2
    echo "Check the available installation with: module spider openblas/0.3.17" >&2
    exit 1
fi

export LD_LIBRARY_PATH="${OPENBLAS_LIB}:${LD_LIBRARY_PATH}"

# ---------------------------------------------------------------------------
# Submit from the HelioCubed repository root:
#
#     sbatch job_convergence_anvil.sh
# ---------------------------------------------------------------------------

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
cd "${PROJECT_DIR}"

export HELIOCUBED_EXE="${PROJECT_DIR}/exec/cubedSphereTest.exe"
export HELIOCUBED_CONVERGENCE_INPUT="${PROJECT_DIR}/exec/inputs_convergence"

RESULTS_ROOT="${PROJECT_DIR}/Convergence_jobs/${SLURM_JOB_ID}"
export HELIOCUBED_CONVERGENCE_RESULTS="${RESULTS_ROOT}"

# ---------------------------------------------------------------------------
# Validate files and shared libraries before launching MPI.
# ---------------------------------------------------------------------------

if [[ ! -x "${HELIOCUBED_EXE}" ]]; then
    echo "ERROR: Executable missing or not executable: ${HELIOCUBED_EXE}" >&2
    exit 1
fi

if [[ ! -r "${HELIOCUBED_CONVERGENCE_INPUT}" ]]; then
    echo "ERROR: Input file missing or unreadable: ${HELIOCUBED_CONVERGENCE_INPUT}" >&2
    exit 1
fi

if [[ ! -r "${PROJECT_DIR}/scripts/run_convergence_test.sh" ]]; then
    echo "ERROR: Convergence driver missing or unreadable." >&2
    exit 1
fi

echo "Checking executable shared libraries:"
LIBRARY_CHECK="$(ldd "${HELIOCUBED_EXE}")"
printf '%s\n' "${LIBRARY_CHECK}"

if [[ "${LIBRARY_CHECK}" == *"not found"* ]]; then
    echo "ERROR: Executable has unresolved shared libraries." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Convergence-test configuration.
# ---------------------------------------------------------------------------

# Problem types:
#   2 = radial pulse
#   3 = radial pulse with constant magnetic field
#   4 = non-radial pulse
#   5 = non-radial pulse with constant magnetic field
export HELIOCUBED_CONVERGENCE_PROBLEM_TYPE=3

# 1 = spatial convergence
# 2 = combined space-and-time convergence
export HELIOCUBED_CONVERGENCE_TEST_TYPE=1

# Base resolution. Levels 1 and 2 use 2x and 4x these values.
export HELIOCUBED_CONVERGENCE_DOMAIN_SIZE=45
export HELIOCUBED_CONVERGENCE_THICKNESS=60

# Fixed box sizes produce 18, 144, and 1152 boxes.
export HELIOCUBED_CONVERGENCE_BOX_SIZE_NONRAD=45
export HELIOCUBED_CONVERGENCE_BOX_SIZE_RAD=20

# Spatial convergence uses this iteration count at every level.
export HELIOCUBED_CONVERGENCE_MAX_ITER=20

# MPI ranks for levels 0, 1, and 2.
# The comparison step uses the level-0 count.
export HELIOCUBED_CONVERGENCE_NPROCS_BY_LEVEL="18 144 1152"

# Time integrator. Spatial convergence uses the same timestep at every level.
export HELIOCUBED_CONVERGENCE_TEMPORAL_ORDER=4

export MPIEXEC=mpirun

# ---------------------------------------------------------------------------
# Diagnostics.
# ---------------------------------------------------------------------------

echo "============================================================"
echo "Slurm job ID:       ${SLURM_JOB_ID}"
echo "Allocated nodes:    ${SLURM_JOB_NUM_NODES}"
echo "Allocated tasks:    ${SLURM_NTASKS}"
echo "Project directory:  ${PROJECT_DIR}"
echo "Executable:         ${HELIOCUBED_EXE}"
echo "Results root:       ${HELIOCUBED_CONVERGENCE_RESULTS}"
echo "MPI ranks:          ${HELIOCUBED_CONVERGENCE_NPROCS_BY_LEVEL}"
echo "============================================================"

# ---------------------------------------------------------------------------
# Run levels 0, 1, and 2 sequentially, followed by level 3 comparison.
# The driver creates the results directories automatically.
# ---------------------------------------------------------------------------

bash "${PROJECT_DIR}/scripts/run_convergence_test.sh"

echo
echo "Convergence job completed successfully."
echo "Results are in:"
echo "  ${HELIOCUBED_CONVERGENCE_RESULTS}/radial_pulse_constant_B"