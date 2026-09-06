#!/usr/bin/env bash

# Run HelioCubed's three-level convergence test for a radial or non-radial
# smooth pulse, with or without constant B. Each resolution is a separate MPI
# launch, allowing the finer levels to use more ranks. A final, inexpensive
# launch reloads the three solutions and reports the observed order for all
# eight variables.
#
# This script takes no command-line settings. Edit USER CONFIGURATION below.
# For automation, the same settings can be overridden with the
# HELIOCUBED_CONVERGENCE_* environment variables shown on each assignment.

set -Eeuo pipefail

# =============================================================================
# USER CONFIGURATION -- edit these values before running the script
# =============================================================================

# Pulse to test:
#   2 = radial pulse
#   3 = radial pulse with constant Cartesian magnetic field
#   4 = non-radial pulse
#   5 = non-radial pulse with constant Cartesian magnetic field
PROBLEM_TYPE="${HELIOCUBED_CONVERGENCE_PROBLEM_TYPE:-5}"

# Convergence mode:
#   1 = spatial convergence (refine the mesh only)
#   2 = space-and-time convergence (refine the mesh and timestep together)
CONVERGENCE_TEST_TYPE="${HELIOCUBED_CONVERGENCE_TEST_TYPE:-1}"

# Base mesh resolution. Both convergence modes also run (2x, 2x) and
# (4x, 4x) versions of these angular and radial resolutions.
DOMAIN_SIZE="${HELIOCUBED_CONVERGENCE_DOMAIN_SIZE:-45}"
THICKNESS="${HELIOCUBED_CONVERGENCE_THICKNESS:-60}"

# Simulation box dimensions. Each must divide its corresponding base mesh
# dimension exactly. They remain fixed as the mesh is refined, creating more
# boxes for additional MPI ranks. The comparison stage redistributes the saved
# solutions onto nested boxes before averaging fine data to coarse grids.
BOX_SIZE_NONRAD="${HELIOCUBED_CONVERGENCE_BOX_SIZE_NONRAD:-45}"
BOX_SIZE_RAD="${HELIOCUBED_CONVERGENCE_BOX_SIZE_RAD:-20}"

# MAX_ITER applies to the base level. Type 1 uses this count on every level;
# type 2 uses 2*MAX_ITER and 4*MAX_ITER on the finer levels.
MAX_ITER="${HELIOCUBED_CONVERGENCE_MAX_ITER:-3}"

# MPI ranks for the base, 2x, and 4x levels. Increase these according to the
# resources in your allocation. The final comparison uses the base count.
if [[ -n "${HELIOCUBED_CONVERGENCE_NPROCS_BY_LEVEL:-}" ]]; then
    read -r -a NPROCS_BY_LEVEL <<< "${HELIOCUBED_CONVERGENCE_NPROCS_BY_LEVEL}"
else
    NPROCS_BY_LEVEL=(18 18 18)
fi

# Supported time integrators are 1, 3, and 4. Fourth order is the normal
# choice when measuring the full space-and-time accuracy of the scheme.
TEMPORAL_ORDER="${HELIOCUBED_CONVERGENCE_TEMPORAL_ORDER:-1}"

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
    2) PROBLEM_NAME="radial_pulse" ;;
    3) PROBLEM_NAME="radial_pulse_constant_B" ;;
    4) PROBLEM_NAME="non_radial_pulse" ;;
    5) PROBLEM_NAME="non_radial_pulse_constant_B" ;;
    *)
        echo "Error: PROBLEM_TYPE must be 2, 3, 4, or 5." >&2
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
if (( ${#NPROCS_BY_LEVEL[@]} != 3 )); then
    echo "Error: NPROCS_BY_LEVEL must contain exactly three process counts." >&2
    exit 1
fi
for level in 0 1 2; do
    require_positive_integer "NPROCS_BY_LEVEL[${level}]" "${NPROCS_BY_LEVEL[${level}]}"
done

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
for level in 0 1 2; do
    scale=$((1 << level))
    level_box_count=$((BASE_BOX_COUNT * scale * scale * scale))
    if (( NPROCS_BY_LEVEL[level] > level_box_count )); then
        echo "Error: NPROCS_BY_LEVEL[${level}] (${NPROCS_BY_LEVEL[${level}]}) exceeds the ${level_box_count} computational boxes at level ${level}." >&2
        echo "Reduce that process count or reduce the box sizes so every MPI rank owns a box." >&2
        exit 1
    fi
    if (( level > 0 )); then
        if (( NPROCS_BY_LEVEL[level] < NPROCS_BY_LEVEL[level - 1] )); then
            echo "Error: NPROCS_BY_LEVEL must be nondecreasing." >&2
            exit 1
        fi
    fi
done
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
    local convergence_level="$2"
    local resolution_scale="$3"
    local convergence_dt="$4"
    local level_domain_size=$((resolution_scale * DOMAIN_SIZE))
    local level_thickness=$((resolution_scale * THICKNESS))
    local level_max_iter="${MAX_ITER}"

    if [[ "${CONVERGENCE_TEST_TYPE}" -eq 2 ]]; then
        level_max_iter=$((resolution_scale * MAX_ITER))
    fi

    awk \
        -v convergence_test_type="${CONVERGENCE_TEST_TYPE}" \
        -v convergence_level="${convergence_level}" \
        -v convergence_dt="${convergence_dt}" \
        -v problem_type="${PROBLEM_TYPE}" \
        -v domain_size="${level_domain_size}" \
        -v thickness="${level_thickness}" \
        -v box_size_nonrad="${BOX_SIZE_NONRAD}" \
        -v box_size_rad="${BOX_SIZE_RAD}" \
        -v max_iter="${level_max_iter}" \
        -v temporal_order="${TEMPORAL_ORDER}" \
        -v data_prefix="${PROBLEM_NAME}_convergence_L${convergence_level}" \
        -v checkpoint_prefix="${PROBLEM_NAME}_convergence_checkpoint_L${convergence_level}" \
        '
        BEGIN { found_level = 0; found_dt = 0 }
        $1 == "-convTestType"           { $2 = convergence_test_type }
        $1 == "-convergence_level"      { $2 = convergence_level; found_level = 1 }
        $1 == "-convergence_dt"         { $2 = sprintf("%.17g", convergence_dt); found_dt = 1 }
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
        END {
            if (!found_level) print "-convergence_level " convergence_level
            if (!found_dt) printf "-convergence_dt %.17g\n", convergence_dt
        }
        ' "${BASE_INPUT}" > "${output_file}"
}

clean_results_directory
LOG_FILE="${RESULTS_DIR}/run.log"
SUMMARY_FILE="${RESULTS_DIR}/convergence_summary.txt"

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
echo "MPI processes by level: ${NPROCS_BY_LEVEL[*]}"
echo "Results: ${RESULTS_DIR}"

BASE_DT=""
for level in 0 1 2; do
    scale=$((1 << level))
    level_domain_size=$((scale * DOMAIN_SIZE))
    level_thickness=$((scale * THICKNESS))
    level_max_iter="${MAX_ITER}"
    level_dt="-1.0"
    if [[ "${CONVERGENCE_TEST_TYPE}" -eq 2 ]]; then
        level_max_iter=$((scale * MAX_ITER))
    fi
    if (( level > 0 )); then
        if [[ "${CONVERGENCE_TEST_TYPE}" -eq 1 ]]; then
            level_dt="${BASE_DT}"
        else
            level_dt="$(awk -v dt="${BASE_DT}" -v divisor="${scale}" 'BEGIN { printf "%.17g", dt / divisor }')"
        fi
    fi

    input_file="${RESULTS_DIR}/inputs_level${level}"
    level_log="${RESULTS_DIR}/level${level}.log"
    make_convergence_input "${input_file}" "${level}" "${scale}" "${level_dt}"

    printf '\n=== Running convergence level %d: %dx%d, %d radial cells, %d MPI processes, %d iterations ===\n' \
        "${level}" "${level_domain_size}" "${level_domain_size}" "${level_thickness}" \
        "${NPROCS_BY_LEVEL[${level}]}" "${level_max_iter}" | tee -a "${LOG_FILE}"
    (
        cd "${RESULTS_DIR}"
        "${MPIEXEC}" -np "${NPROCS_BY_LEVEL[${level}]}" "${EXECUTABLE}" "${input_file}"
    ) 2>&1 | tee -a "${LOG_FILE}" "${level_log}"

    if (( level == 0 )); then
        BASE_DT="$(awk '$1 == "Convergence" && $2 == "timestep" && $3 == "=" { value = $4 } END { print value }' "${level_log}")"
        if ! [[ "${BASE_DT}" =~ ^[+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]] \
            || ! awk -v dt="${BASE_DT}" 'BEGIN { exit !(dt > 0) }'; then
            echo "Error: could not read a positive base convergence timestep from ${level_log}." >&2
            exit 1
        fi
        echo "Base convergence timestep: ${BASE_DT}" | tee -a "${LOG_FILE}"
    fi
done

COMPARE_INPUT="${RESULTS_DIR}/inputs_compare"
COMPARE_LOG="${RESULTS_DIR}/compare.log"
make_convergence_input "${COMPARE_INPUT}" 3 1 -1.0
printf '\n=== Comparing the three convergence levels with %d MPI processes ===\n' \
    "${NPROCS_BY_LEVEL[0]}" | tee -a "${LOG_FILE}"
(
    cd "${RESULTS_DIR}"
    "${MPIEXEC}" -np "${NPROCS_BY_LEVEL[0]}" "${EXECUTABLE}" "${COMPARE_INPUT}"
) 2>&1 | tee -a "${LOG_FILE}" "${COMPARE_LOG}"

awk '/Lev = [01], component = [0-7], error =/ || /order of accuracy for var [0-7] =/ { print }' \
    "${COMPARE_LOG}" > "${SUMMARY_FILE}"

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
