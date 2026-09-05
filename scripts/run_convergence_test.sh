#!/usr/bin/env bash

# Run HelioCubed's built-in three-level convergence test for either the radial
# or non-radial smooth pulse. The solver compares its base, 2x, and 4x
# solutions and reports the observed order for all eight variables.
#
# This script takes no command-line settings. Edit USER CONFIGURATION below.

set -Eeuo pipefail

# =============================================================================
# USER CONFIGURATION -- edit these values before running the script
# =============================================================================

# Pulse to test:
#   0 = smooth radial pulse
#   1 = smooth non-radial pulse
PROBLEM_TYPE=1

# Convergence mode:
#   1 = spatial convergence (refine the mesh only)
#   2 = space-and-time convergence (refine the mesh and timestep together)
CONVERGENCE_TEST_TYPE=1

# Base mesh resolution. Both convergence modes also run (2x, 2x) and
# (4x, 4x) versions of these angular and radial resolutions.
DOMAIN_SIZE=60
THICKNESS=90

# Base computational/MPI box dimensions. Each must divide its corresponding
# base mesh dimension exactly. The solver scales the boxes with each level.
BOX_SIZE_NONRAD=60
BOX_SIZE_RAD=30

# MAX_ITER applies to the base level. Type 1 uses this count on every level;
# type 2 uses 2*MAX_ITER and 4*MAX_ITER on the finer levels.
MAX_ITER=3
NPROCS=18

# Supported time integrators are 1, 3, and 4. Fourth order is the normal
# choice when measuring the full space-and-time accuracy of the scheme.
TEMPORAL_ORDER=1

# =============================================================================
# END USER CONFIGURATION
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
EXECUTABLE="${HELIOCUBED_EXE:-${ROOT_DIR}/exec/cubedSphereTest.exe}"
BASE_INPUT="${HELIOCUBED_CONVERGENCE_INPUT:-${ROOT_DIR}/exec/inputs_convergence}"
RESULTS_ROOT="${HELIOCUBED_CONVERGENCE_RESULTS:-${ROOT_DIR}/Convergence_results}"
MPIEXEC="${MPIEXEC:-mpirun}"

require_positive_integer() {
    local variable_name="$1"
    local value="$2"

    if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: ${variable_name} must be a positive integer (got '${value}')." >&2
        exit 1
    fi
}

case "${PROBLEM_TYPE}" in
    0) PROBLEM_NAME="radial_pulse" ;;
    1) PROBLEM_NAME="non_radial_pulse" ;;
    *)
        echo "Error: PROBLEM_TYPE must be 0 (radial) or 1 (non-radial)." >&2
        exit 1
        ;;
esac

case "${CONVERGENCE_TEST_TYPE}" in
    1) CONVERGENCE_MODE="space only" ;;
    2) CONVERGENCE_MODE="space and time" ;;
    *)
        echo "Error: CONVERGENCE_TEST_TYPE must be 1 (space) or 2 (space and time)." >&2
        exit 1
        ;;
esac

require_positive_integer "DOMAIN_SIZE" "${DOMAIN_SIZE}"
require_positive_integer "THICKNESS" "${THICKNESS}"
require_positive_integer "BOX_SIZE_NONRAD" "${BOX_SIZE_NONRAD}"
require_positive_integer "BOX_SIZE_RAD" "${BOX_SIZE_RAD}"
require_positive_integer "MAX_ITER" "${MAX_ITER}"
require_positive_integer "NPROCS" "${NPROCS}"

if (( DOMAIN_SIZE % BOX_SIZE_NONRAD != 0 )); then
    echo "Error: BOX_SIZE_NONRAD (${BOX_SIZE_NONRAD}) must divide DOMAIN_SIZE (${DOMAIN_SIZE})." >&2
    exit 1
fi
if (( THICKNESS % BOX_SIZE_RAD != 0 )); then
    echo "Error: BOX_SIZE_RAD (${BOX_SIZE_RAD}) must divide THICKNESS (${THICKNESS})." >&2
    exit 1
fi

BASE_BOX_COUNT=$((
    6 * (DOMAIN_SIZE / BOX_SIZE_NONRAD) * (DOMAIN_SIZE / BOX_SIZE_NONRAD)
    * (THICKNESS / BOX_SIZE_RAD)
))
if (( NPROCS > BASE_BOX_COUNT )); then
    echo "Error: NPROCS (${NPROCS}) exceeds the ${BASE_BOX_COUNT} computational boxes." >&2
    echo "Reduce NPROCS or reduce the box sizes so every MPI rank owns a box." >&2
    exit 1
fi
case "${TEMPORAL_ORDER}" in
    1|3|4) ;;
    *)
        echo "Error: TEMPORAL_ORDER must be 1, 3, or 4 (got '${TEMPORAL_ORDER}')." >&2
        exit 1
        ;;
esac

if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "Error: executable not found or not executable: ${EXECUTABLE}" >&2
    echo "Build it with: make -C \"${ROOT_DIR}/exec\" -f GNUmakefile_TS -j4" >&2
    exit 1
fi
if [[ ! -f "${BASE_INPUT}" ]]; then
    echo "Error: convergence input template not found: ${BASE_INPUT}" >&2
    exit 1
fi
if ! command -v "${MPIEXEC}" >/dev/null 2>&1; then
    echo "Error: MPI launcher not found: ${MPIEXEC}" >&2
    exit 1
fi

mkdir -p "${RESULTS_ROOT}"
RESULTS_ROOT="$(cd "${RESULTS_ROOT}" && pwd -P)"
RESULTS_DIR="${RESULTS_ROOT}/${PROBLEM_NAME}"

clean_results_directory() {
    local expected_results_dir="${RESULTS_ROOT}/${PROBLEM_NAME}"

    if [[ -z "${RESULTS_DIR}" || "${RESULTS_DIR}" == "/" \
        || "${RESULTS_DIR}" == "${RESULTS_ROOT}" \
        || "${RESULTS_DIR}" != "${expected_results_dir}" ]]; then
        echo "Error: refusing to clean unsafe results directory: '${RESULTS_DIR}'." >&2
        exit 1
    fi

    if [[ -e "${RESULTS_DIR}" || -L "${RESULTS_DIR}" ]]; then
        echo "Cleaning previous convergence results: ${RESULTS_DIR}"
        rm -rf -- "${RESULTS_DIR}"
    fi
    mkdir -p "${RESULTS_DIR}"
}

make_convergence_input() {
    local output_file="$1"

    awk \
        -v convergence_test_type="${CONVERGENCE_TEST_TYPE}" \
        -v problem_type="${PROBLEM_TYPE}" \
        -v domain_size="${DOMAIN_SIZE}" \
        -v thickness="${THICKNESS}" \
        -v box_size_nonrad="${BOX_SIZE_NONRAD}" \
        -v box_size_rad="${BOX_SIZE_RAD}" \
        -v max_iter="${MAX_ITER}" \
        -v temporal_order="${TEMPORAL_ORDER}" \
        -v data_prefix="${PROBLEM_NAME}_convergence" \
        -v checkpoint_prefix="${PROBLEM_NAME}_convergence_checkpoint" \
        '
        $1 == "-convTestType"           { $2 = convergence_test_type }
        $1 == "-init_condition_type"    { $2 = problem_type }
        $1 == "-domainSize"             { $2 = domain_size }
        $1 == "-thickness"              { $2 = thickness }
        $1 == "-boxSize_nonrad"         { $2 = box_size_nonrad }
        $1 == "-boxSize_rad"            { $2 = box_size_rad }
        $1 == "-max_iter"               { $2 = max_iter }
        $1 == "-write_cadence"          { $2 = 1000000000 }
        $1 == "-checkpoint_cadence"     { $2 = 1000000000 }
        $1 == "-P_floor_cadence"        { $2 = 1000000000 }
        $1 == "-data_file_prefix"       { $2 = data_prefix }
        $1 == "-checkpoint_file_prefix" { $2 = checkpoint_prefix }
        $1 == "-temporal_order"         { $2 = temporal_order }
        $1 == "-radial_refinement"      { $2 = 0 }
        { print }
        ' "${BASE_INPUT}" > "${output_file}"
}

clean_results_directory
INPUT_FILE="${RESULTS_DIR}/inputs"
LOG_FILE="${RESULTS_DIR}/run.log"
SUMMARY_FILE="${RESULTS_DIR}/convergence_summary.txt"
make_convergence_input "${INPUT_FILE}"

FINE_DOMAIN_SIZE=$((4 * DOMAIN_SIZE))
FINE_THICKNESS=$((4 * THICKNESS))

echo "Problem: ${PROBLEM_NAME} (init_condition_type=${PROBLEM_TYPE})"
echo "Convergence mode: convTestType=${CONVERGENCE_TEST_TYPE} (${CONVERGENCE_MODE})"
echo "Angular resolutions: ${DOMAIN_SIZE}, $((2 * DOMAIN_SIZE)), ${FINE_DOMAIN_SIZE}"
echo "Radial resolutions: ${THICKNESS}, $((2 * THICKNESS)), ${FINE_THICKNESS}"
if [[ "${CONVERGENCE_TEST_TYPE}" -eq 2 ]]; then
    echo "Iteration counts: ${MAX_ITER}, $((2 * MAX_ITER)), $((4 * MAX_ITER))"
else
    echo "Iteration counts: ${MAX_ITER}, ${MAX_ITER}, ${MAX_ITER}"
fi
echo "MPI processes: ${NPROCS}"
echo "Results: ${RESULTS_DIR}"

(
    cd "${RESULTS_DIR}"
    "${MPIEXEC}" -np "${NPROCS}" "${EXECUTABLE}" "${INPUT_FILE}"
) 2>&1 | tee "${LOG_FILE}"

awk '/Lev = [01], component = [0-7], error =/ || /order of accuracy for var [0-7] =/ { print }' \
    "${LOG_FILE}" > "${SUMMARY_FILE}"

if [[ ! -s "${SUMMARY_FILE}" ]]; then
    echo "Error: the solver completed without reporting convergence results." >&2
    exit 1
fi

echo
echo "=== Convergence summary ==="
cat "${SUMMARY_FILE}"
echo
echo "Component map: 0=density, 1-3=momentum, 4=energy, 5-7=magnetic field"
echo "Full log: ${LOG_FILE}"
echo "Summary: ${SUMMARY_FILE}"
