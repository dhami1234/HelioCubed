#!/bin/bash
#SBATCH --account=mca07s033
#SBATCH --partition=wholenode
#SBATCH --nodes=9
#SBATCH --ntasks=1152
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --job-name=hc-conv
#SBATCH --output=hc-conv.o%j
#SBATCH --error=hc-conv.e%j

# Submit one independent Slurm job per selected convergence problem:
#
#   bash scripts/run_convergence_test_anvil.sh
#   HELIOCUBED_CONVERGENCE_TEST_CASES="2 3 4 5" bash scripts/run_convergence_test_anvil.sh
#
# Run on an Anvil login node. Each job runs levels 0, 1, 2, then comparison.
# Slurm can run the problem jobs concurrently as resources become available.
# This script submits itself in worker mode; do not invoke --run-case manually.

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: line ${LINENO}: ${BASH_COMMAND} (exit ${rc})" >&2; exit "$rc"' ERR

run_anvil_case() {
    local project_dir="$1"
    export HELIOCUBED_CONVERGENCE_PROBLEM_TYPE="$2"
    local max_time="$3"
    local results_root="${HELIOCUBED_CONVERGENCE_RESULTS}/${SLURM_JOB_ID}"

    # Compiler must be loaded before MPI.
    module purge
    module load xalt/2.10.45
    module load intel/19.0.5.281
    module load numactl/2.0.14
    module load libszip/2.1.1
    module load zlib/1.2.11
    module load mvapich2/2.3.6
    module load hdf5/1.10.7

    module list
    command -v "${MPIEXEC}"

    export MV2_HOMOGENEOUS_CLUSTER=1
    export OMP_NUM_THREADS=1
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close

    # ---------------------------------------------------------------------------
    # Shared-library paths.
    # ---------------------------------------------------------------------------

    export LD_LIBRARY_PATH="${HELIOCUBED_ANVIL_HDF5_LIB:-/apps/spack/anvil/apps/hdf5/1.10.7-intel-19.0.5-gddu73a/lib}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
    export LD_LIBRARY_PATH="${HELIOCUBED_ANVIL_GSL_LIB:-/home/x-dhami/gsl/lib}:${LD_LIBRARY_PATH}"

    # OpenBLAS location specified in exec/GNUmakefile_Anvil.
    local OPENBLAS_LIB="${HELIOCUBED_ANVIL_OPENBLAS_LIB:-/apps/spack/anvil/apps/openblas/0.3.17-gcc-11.2.0-2qrsari/lib}"

    if [[ ! -r "${OPENBLAS_LIB}/libopenblas.so.0" ]]; then
        echo "ERROR: Cannot find ${OPENBLAS_LIB}/libopenblas.so.0" >&2
        echo "Check the available installation with: module spider openblas/0.3.17" >&2
        exit 1
    fi

    export LD_LIBRARY_PATH="${OPENBLAS_LIB}:${LD_LIBRARY_PATH}"

    cd "${project_dir}"
    echo "Checking executable shared libraries:"
    local library_check
    library_check="$(ldd "${HELIOCUBED_EXE}")"
    printf '%s\n' "${library_check}"
    if [[ "${library_check}" == *"not found"* ]]; then
        echo "ERROR: Executable has unresolved shared libraries." >&2
        exit 1
    fi

    # A job-specific template keeps max_time identical at every resolution
    # without modifying the shared input file or another problem's settings.
    mkdir -p "${results_root}"
    local case_input="${results_root}/inputs_convergence"
    awk -v max_time="${max_time}" '
        $1 == "-max_time" { $2 = max_time; found = 1 }
        { print }
        END { if (!found) print "-max_time " max_time }
    ' "${HELIOCUBED_CONVERGENCE_INPUT}" > "${case_input}"
    export HELIOCUBED_CONVERGENCE_INPUT="${case_input}"
    export HELIOCUBED_CONVERGENCE_RESULTS="${results_root}"

    echo "Slurm job ID:      ${SLURM_JOB_ID}"
    echo "Problem type:      ${HELIOCUBED_CONVERGENCE_PROBLEM_TYPE}"
    echo "Maximum time:      ${max_time}"
    echo "Allocated tasks:   ${SLURM_NTASKS}"
    echo "MPI ranks:         ${HELIOCUBED_CONVERGENCE_NPROCS_BY_LEVEL}"
    echo "Results root:      ${results_root}"

    bash "${project_dir}/scripts/run_convergence_test.sh"
    echo "Convergence job completed successfully. Results: ${results_root}"
}

if [[ "${1:-}" == "--run-case" ]]; then
    if [[ "$#" -ne 4 ]]; then
        echo "ERROR: internal worker invocation requires a project directory, problem ID, and max_time." >&2
        exit 2
    fi
    run_anvil_case "$2" "$3" "$4"
    exit 0
fi
if [[ "$#" -ne 0 ]]; then
    echo "ERROR: edit USER CONFIGURATION or use the documented environment overrides." >&2
    exit 2
fi

# =============================================================================
# USER CONFIGURATION -- edit these values before submitting
# =============================================================================

# Select one or more problems, e.g. TEST_CASES=(2 3 4 5).
#   0 = file-driven solar wind (use an absolute -BC_file path)
#   1 = analytic outflow
#   2 = radial pulse
#   3 = radial pulse with constant Cartesian magnetic field
#   4 = non-radial pulse
#   5 = non-radial pulse with constant Cartesian magnetic field
#  10 = smooth exact expanding wind, B=0
#  11 = the same wind with uniform Cartesian B(t)
# Both require wind_boundary_order 4 and Sun_gravity 0.
TEST_CASES=(3 4 5)
if [[ -n "${HELIOCUBED_CONVERGENCE_TEST_CASES:-}" ]]; then
    read -r -a TEST_CASES <<< "${HELIOCUBED_CONVERGENCE_TEST_CASES}"
elif [[ -n "${HELIOCUBED_CONVERGENCE_PROBLEM_TYPE:-}" ]]; then
    # Preserve the single-problem environment setting.
    TEST_CASES=("${HELIOCUBED_CONVERGENCE_PROBLEM_TYPE}")
fi

# Maximum simulation time in seconds, independently for each problem.
CASE_0_MAX_TIME="${HELIOCUBED_CONVERGENCE_CASE_0_MAX_TIME:-36000.0}"
CASE_1_MAX_TIME="${HELIOCUBED_CONVERGENCE_CASE_1_MAX_TIME:-36000.0}"
# Two initial cycles become one at tau=29928 s for the standard inner radius.
CASE_10_MAX_TIME="${HELIOCUBED_CONVERGENCE_CASE_10_MAX_TIME:-29928.0}"
CASE_11_MAX_TIME="${HELIOCUBED_CONVERGENCE_CASE_11_MAX_TIME:-29928.0}"
CASE_2_MAX_TIME="${HELIOCUBED_CONVERGENCE_CASE_2_MAX_TIME:-47600.0}"
CASE_3_MAX_TIME="${HELIOCUBED_CONVERGENCE_CASE_3_MAX_TIME:-47600.0}"
CASE_4_MAX_TIME="${HELIOCUBED_CONVERGENCE_CASE_4_MAX_TIME:-360000.0}"
CASE_5_MAX_TIME="${HELIOCUBED_CONVERGENCE_CASE_5_MAX_TIME:-36000.0}"

# 1 = spatial convergence; 2 = combined space-and-time convergence.
export HELIOCUBED_CONVERGENCE_TEST_TYPE="${HELIOCUBED_CONVERGENCE_TEST_TYPE:-1}"

# Base mesh. Levels 1 and 2 use 2x and 4x these dimensions.
export HELIOCUBED_CONVERGENCE_DOMAIN_SIZE="${HELIOCUBED_CONVERGENCE_DOMAIN_SIZE:-45}"
export HELIOCUBED_CONVERGENCE_THICKNESS="${HELIOCUBED_CONVERGENCE_THICKNESS:-60}"
# Fixed box sizes produce 18, 144, and 1152 boxes at the default mesh.
export HELIOCUBED_CONVERGENCE_BOX_SIZE_NONRAD="${HELIOCUBED_CONVERGENCE_BOX_SIZE_NONRAD:-45}"
export HELIOCUBED_CONVERGENCE_BOX_SIZE_RAD="${HELIOCUBED_CONVERGENCE_BOX_SIZE_RAD:-20}"

# The earlier of max_iter and max_time wins. Increase this count if you want
# runs to reach max_time. Mode 1 uses this count at every level; mode 2 uses
# this count, 2x this count, and 4x this count.
export HELIOCUBED_CONVERGENCE_MAX_ITER="${HELIOCUBED_CONVERGENCE_MAX_ITER:-200000}"
# MPI ranks for levels 0, 1, 2. Comparison uses the level-0 count.
export HELIOCUBED_CONVERGENCE_NPROCS_BY_LEVEL="${HELIOCUBED_CONVERGENCE_NPROCS_BY_LEVEL:-18 144 1152}"
# Leave empty to preserve -temporal_order in the input template.
export HELIOCUBED_CONVERGENCE_TEMPORAL_ORDER="${HELIOCUBED_CONVERGENCE_TEMPORAL_ORDER:-}"
export MPIEXEC="${MPIEXEC:-mpirun}"

# Slurm resources PER PROBLEM, not shared across the selected problems.
NODES="${HELIOCUBED_ANVIL_NODES:-9}"
NPROCS="${HELIOCUBED_ANVIL_NTASKS:-1152}"
WALLTIME="${HELIOCUBED_ANVIL_WALLTIME:-72:00:00}"
PARTITION="${HELIOCUBED_ANVIL_PARTITION:-wholenode}"
ACCOUNT="${HELIOCUBED_ANVIL_ACCOUNT:-mca07s033}"

# =============================================================================
# END USER CONFIGURATION
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
export HELIOCUBED_EXE="${HELIOCUBED_EXE:-${PROJECT_DIR}/exec/cubedSphereTest.exe}"
export HELIOCUBED_CONVERGENCE_INPUT="${HELIOCUBED_CONVERGENCE_INPUT:-${PROJECT_DIR}/exec/inputs_convergence}"
# Each worker adds its unique Slurm job ID, preventing simultaneous runs
# (including repeat submissions of the same problem) from sharing outputs.
export HELIOCUBED_CONVERGENCE_RESULTS="${HELIOCUBED_CONVERGENCE_RESULTS:-${PROJECT_DIR}/Convergence_jobs}"
# Resolve relative overrides against the launch directory before Slurm chdir.
for variable in HELIOCUBED_EXE HELIOCUBED_CONVERGENCE_INPUT HELIOCUBED_CONVERGENCE_RESULTS; do
    if [[ "${!variable}" != /* ]]; then
        printf -v "${variable}" '%s/%s' "${PWD}" "${!variable}"
    fi
done

require_positive_integer() {
    if ! [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: $1 must be a positive integer (got '$2')." >&2
        exit 1
    fi
}

if (( ${#TEST_CASES[@]} == 0 )); then
    echo "ERROR: TEST_CASES must contain at least one problem ID." >&2
    exit 1
fi
seen_case_ids=" "
for case_id in "${TEST_CASES[@]}"; do
    case "${case_id}" in
        0|1|2|3|4|5|10|11) ;;
        *) echo "ERROR: unsupported problem ID '${case_id}'; choose 0, 1, 2, 3, 4, 5, 10, or 11." >&2; exit 1 ;;
    esac
    if [[ "${seen_case_ids}" == *" ${case_id} "* ]]; then
        echo "ERROR: duplicate problem ID '${case_id}' in TEST_CASES." >&2
        exit 1
    fi
    seen_case_ids+="${case_id} "
    time_variable="CASE_${case_id}_MAX_TIME"
    if ! awk -v value="${!time_variable}" 'BEGIN {
        exit !(value ~ /^[+]?[0-9]*[.]?[0-9]+([eE][+-]?[0-9]+)?$/ && value+0 > 0)
    }'; then
        echo "ERROR: ${time_variable} must be a positive number (got '${!time_variable}')." >&2
        exit 1
    fi
done
require_positive_integer NODES "${NODES}"
require_positive_integer NPROCS "${NPROCS}"
read -r -a ranks_by_level <<< "${HELIOCUBED_CONVERGENCE_NPROCS_BY_LEVEL}"
if (( ${#ranks_by_level[@]} != 3 )); then
    echo "ERROR: HELIOCUBED_CONVERGENCE_NPROCS_BY_LEVEL must contain three counts." >&2
    exit 1
fi
for ranks in "${ranks_by_level[@]}"; do
    require_positive_integer MPI_RANKS "${ranks}"
    if (( ranks > NPROCS )); then
        echo "ERROR: requested ${ranks} MPI ranks exceeds the per-job allocation (${NPROCS})." >&2
        exit 1
    fi
done
if ! [[ "${WALLTIME}" =~ ^[0-9]+:[0-5][0-9]:[0-5][0-9]$ ]]; then
    echo "ERROR: WALLTIME must use hh:mm:ss format." >&2
    exit 1
fi
if [[ ! -x "${HELIOCUBED_EXE}" || ! -r "${HELIOCUBED_CONVERGENCE_INPUT}" \
    || ! -r "${PROJECT_DIR}/scripts/run_convergence_test.sh" ]]; then
    echo "ERROR: executable, convergence input, or convergence driver missing/unreadable." >&2
    exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is not available. Run this launcher on an Anvil login node." >&2
    exit 1
fi

mkdir -p "${HELIOCUBED_CONVERGENCE_RESULTS}"
echo "Selected problems: ${TEST_CASES[*]}"
echo "Per problem: nodes=${NODES}, tasks=${NPROCS}, time=${WALLTIME}"
echo "Results: ${HELIOCUBED_CONVERGENCE_RESULTS}/<job_id>/<problem_name>"
for case_id in "${TEST_CASES[@]}"; do
    time_variable="CASE_${case_id}_MAX_TIME"
    submission="$(sbatch \
        --export=ALL \
        --account="${ACCOUNT}" \
        --partition="${PARTITION}" \
        --nodes="${NODES}" \
        --ntasks="${NPROCS}" \
        --time="${WALLTIME}" \
        --job-name="hc-conv-${case_id}" \
        --output="${HELIOCUBED_CONVERGENCE_RESULTS}/hc-conv-${case_id}.o%j" \
        --error="${HELIOCUBED_CONVERGENCE_RESULTS}/hc-conv-${case_id}.e%j" \
        --chdir="${PROJECT_DIR}" \
        "${SCRIPT_PATH}" --run-case "${PROJECT_DIR}" "${case_id}" "${!time_variable}")"
    echo "Problem ${case_id}, max_time=${!time_variable}: ${submission}"
done
echo "Submitted ${#TEST_CASES[@]} convergence job(s)."
echo "Use 'squeue -u ${USER:-$(id -un)}' to monitor them."
