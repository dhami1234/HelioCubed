#!/usr/bin/env bash

# Run selected idealized HelioCubed test problems and plot their z=0 slices.
#
# Usage:
#   ./scripts/run_test_problems.sh [results_directory]
#
# Defaults (set in USER CONFIGURATION below):
#   - Runs the selected idealized test problem(s).
#   - Uses 18 MPI processes, 1000 iterations, and writes z=0 slices every
#     50 iterations (including the initial and final states).
#   - Writes each case into its own directory under ./Test_results.
#   - Deletes an existing result directory for each selected case before
#     starting it, so stale files cannot be mixed into the new run.
#   - Plots density, velocity, pressure, magnetic field, and derived
#     temperature in raw solver/CGS units under each case's plots_z0_raw/.
#
# Resolution:
#   The approximate cell count is:
#
#     6 * domainSize * domainSize * thickness
#
#   boxSize_nonrad and boxSize_rad change domain decomposition, not physical
#   resolution. Their values should divide domainSize and thickness cleanly.
#   NPROCS changes parallelism only; it does not change the mesh resolution.
#
# Useful overrides:
#   PYTHON_BIN=/path/to/python ./scripts/run_test_problems.sh
#   HELIOCUBED_EXE=/path/to/exe ./scripts/run_test_problems.sh /path/to/results
#   HELIOCUBED_TEST_INPUT=/path/to/inputs ./scripts/run_test_problems.sh
#   HELIOCUBED_TEST_CASES="2 7" HELIOCUBED_TEST_MAX_ITER=10 ./scripts/run_test_problems.sh
#   All USER CONFIGURATION settings accept HELIOCUBED_TEST_* overrides below.
#
# The selected Python must provide NumPy and Matplotlib. For every selected
# test, its entire old result directory (including data, logs, and plots) is
# permanently removed immediately before the new run. Unselected test folders
# are left untouched. The script stops if a simulation or plotting step fails.

set -Eeuo pipefail

# =============================================================================
# USER CONFIGURATION -- edit these values before running the script
# =============================================================================

# Tests to run. Use one or more IDs, for example TEST_CASES=(4) or (2 3).
#   2 = radial pulse
#   3 = radial pulse with constant Cartesian magnetic field
#   4 = non-radial pulse
#   5 = non-radial pulse with constant Cartesian magnetic field
#   6 = strong spherical shock
#   7 = strong spherical shock with constant Cartesian magnetic field
#   8 = off-center hydrodynamic blast
#   9 = off-center blast with constant Cartesian B = (0.1, 0, 0) G
TEST_CASES=(8 9 3 5)
if [[ -n "${HELIOCUBED_TEST_CASES:-}" ]]; then
    read -r -a TEST_CASES <<< "${HELIOCUBED_TEST_CASES}"
fi

# Mesh resolution. DOMAIN_SIZE is the angular resolution on each face;
# THICKNESS is the radial resolution.
DOMAIN_SIZE="${HELIOCUBED_TEST_DOMAIN_SIZE:-60}"
THICKNESS="${HELIOCUBED_TEST_THICKNESS:-120}"

# Computational/MPI box dimensions. These must divide the corresponding mesh
# dimensions exactly.
BOX_SIZE_NONRAD="${HELIOCUBED_TEST_BOX_SIZE_NONRAD:-60}"
BOX_SIZE_RAD="${HELIOCUBED_TEST_BOX_SIZE_RAD:-40}"

# Runtime and output settings.
NPROCS="${HELIOCUBED_TEST_NPROCS:-18}"
MAX_ITER="${HELIOCUBED_TEST_MAX_ITER:-100}"
# Maximum simulation time for each test problem. Cases 8 and 9 use the FLASH
# t=0.01 comparison time, converted for the default shell dimensions. Adjust
# their times if rebuilding the shell. HELIOCUBED_TEST_BLAST_MAX_TIME remains
# a shared fallback for both blast cases for backward compatibility.
CASE_2_MAX_TIME="${HELIOCUBED_TEST_CASE_2_MAX_TIME:-3600000.0}" # radial pulse
CASE_3_MAX_TIME="${HELIOCUBED_TEST_CASE_3_MAX_TIME:-3600000.0}" # radial pulse, constant B
CASE_4_MAX_TIME="${HELIOCUBED_TEST_CASE_4_MAX_TIME:-3600000.0}" # non-radial pulse
CASE_5_MAX_TIME="${HELIOCUBED_TEST_CASE_5_MAX_TIME:-3600000.0}" # non-radial pulse, constant B
CASE_6_MAX_TIME="${HELIOCUBED_TEST_CASE_6_MAX_TIME:-3600000.0}" # strong spherical shock
CASE_7_MAX_TIME="${HELIOCUBED_TEST_CASE_7_MAX_TIME:-3600000.0}" # strong shock, constant B
CASE_8_MAX_TIME="${HELIOCUBED_TEST_CASE_8_MAX_TIME:-${HELIOCUBED_TEST_BLAST_MAX_TIME:-5001.23}}" # hydro blast
CASE_9_MAX_TIME="${HELIOCUBED_TEST_CASE_9_MAX_TIME:-${HELIOCUBED_TEST_BLAST_MAX_TIME:-5001.23}}" # MHD blast
BLAST_MAX_ITER="${HELIOCUBED_TEST_MAX_ITER:-10000}"
SLICE_CADENCE="${HELIOCUBED_TEST_SLICE_CADENCE:-50}"
PLOT_WORKERS="${HELIOCUBED_TEST_PLOT_WORKERS:-8}"

# =============================================================================
# END USER CONFIGURATION
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
EXECUTABLE="${HELIOCUBED_EXE:-${ROOT_DIR}/exec/cubedSphereTest.exe}"
BASE_INPUT="${HELIOCUBED_TEST_INPUT:-${ROOT_DIR}/exec/inputs_Test_Problems}"
RESULTS_DIR="${1:-${ROOT_DIR}/Test_results}"
MPIEXEC="${MPIEXEC:-mpirun}"
PLOTTER="${SCRIPT_DIR}/Plot_Test_Slices.py"

case_name_for_id() {
    case "$1" in
        2) printf '%s\n' "02_radial_pulse" ;;
        3) printf '%s\n' "03_radial_pulse_constant_B" ;;
        4) printf '%s\n' "04_non_radial_pulse" ;;
        5) printf '%s\n' "05_non_radial_pulse_constant_B" ;;
        6) printf '%s\n' "06_strong_spherical_shock" ;;
        7) printf '%s\n' "07_strong_spherical_shock_constant_B" ;;
        8) printf '%s\n' "08_hydro_blast" ;;
        9) printf '%s\n' "09_mhd_blast_constant_B" ;;
        *) return 1 ;;
    esac
}

max_time_for_id() {
    case "$1" in
        2) printf '%s\n' "${CASE_2_MAX_TIME}" ;;
        3) printf '%s\n' "${CASE_3_MAX_TIME}" ;;
        4) printf '%s\n' "${CASE_4_MAX_TIME}" ;;
        5) printf '%s\n' "${CASE_5_MAX_TIME}" ;;
        6) printf '%s\n' "${CASE_6_MAX_TIME}" ;;
        7) printf '%s\n' "${CASE_7_MAX_TIME}" ;;
        8) printf '%s\n' "${CASE_8_MAX_TIME}" ;;
        9) printf '%s\n' "${CASE_9_MAX_TIME}" ;;
        *) return 1 ;;
    esac
}

require_positive_integer() {
    local variable_name="$1"
    local value="$2"

    if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: ${variable_name} must be a positive integer (got '${value}')." >&2
        exit 1
    fi
}

require_positive_number() {
    local variable_name="$1"
    local value="$2"

    if ! awk -v value="${value}" 'BEGIN {
        exit !(value ~ /^[0-9]+([.][0-9]*)?([eE][+-]?[0-9]+)?$/ && value+0 > 0)
    }'; then
        echo "Error: ${variable_name} must be a positive number (got '${value}')." >&2
        exit 1
    fi
}

if [[ "${#TEST_CASES[@]}" -eq 0 ]]; then
    echo "Error: TEST_CASES must contain at least one test ID." >&2
    exit 1
fi

CASE_IDS=("${TEST_CASES[@]}")
CASE_NAMES=()
for case_id in "${CASE_IDS[@]}"; do
    if ! case_name="$(case_name_for_id "${case_id}")"; then
        echo "Error: unsupported test ID '${case_id}' in TEST_CASES; choose from 2 through 9." >&2
        exit 1
    fi
    CASE_NAMES+=("${case_name}")
done
TOTAL_CASES="${#CASE_IDS[@]}"

require_positive_integer "DOMAIN_SIZE" "${DOMAIN_SIZE}"
require_positive_integer "THICKNESS" "${THICKNESS}"
require_positive_integer "BOX_SIZE_NONRAD" "${BOX_SIZE_NONRAD}"
require_positive_integer "BOX_SIZE_RAD" "${BOX_SIZE_RAD}"
require_positive_integer "NPROCS" "${NPROCS}"
require_positive_integer "MAX_ITER" "${MAX_ITER}"
require_positive_integer "BLAST_MAX_ITER" "${BLAST_MAX_ITER}"
require_positive_number "CASE_2_MAX_TIME" "${CASE_2_MAX_TIME}"
require_positive_number "CASE_3_MAX_TIME" "${CASE_3_MAX_TIME}"
require_positive_number "CASE_4_MAX_TIME" "${CASE_4_MAX_TIME}"
require_positive_number "CASE_5_MAX_TIME" "${CASE_5_MAX_TIME}"
require_positive_number "CASE_6_MAX_TIME" "${CASE_6_MAX_TIME}"
require_positive_number "CASE_7_MAX_TIME" "${CASE_7_MAX_TIME}"
require_positive_number "CASE_8_MAX_TIME" "${CASE_8_MAX_TIME}"
require_positive_number "CASE_9_MAX_TIME" "${CASE_9_MAX_TIME}"
require_positive_integer "SLICE_CADENCE" "${SLICE_CADENCE}"
require_positive_integer "PLOT_WORKERS" "${PLOT_WORKERS}"

if (( DOMAIN_SIZE % BOX_SIZE_NONRAD != 0 )); then
    echo "Error: BOX_SIZE_NONRAD (${BOX_SIZE_NONRAD}) must divide DOMAIN_SIZE (${DOMAIN_SIZE})." >&2
    exit 1
fi
if (( THICKNESS % BOX_SIZE_RAD != 0 )); then
    echo "Error: BOX_SIZE_RAD (${BOX_SIZE_RAD}) must divide THICKNESS (${THICKNESS})." >&2
    exit 1
fi

if [[ -z "${PYTHON_BIN:-}" ]]; then
    PYTHON_BIN=""
    python_candidates=(
        "${ROOT_DIR}/.venv/bin/python"
        python3
        python
        /opt/homebrew/bin/python3
        /usr/local/bin/python3
    )
    if [[ -n "${CONDA_PREFIX:-}" ]]; then
        python_candidates=("${CONDA_PREFIX}/bin/python" "${python_candidates[@]}")
    fi
    if [[ -n "${VIRTUAL_ENV:-}" ]]; then
        python_candidates=("${VIRTUAL_ENV}/bin/python" "${python_candidates[@]}")
    fi

    for python_candidate in "${python_candidates[@]}"; do
        if command -v "${python_candidate}" >/dev/null 2>&1 \
            && "${python_candidate}" -c 'import matplotlib, numpy' >/dev/null 2>&1; then
            PYTHON_BIN="${python_candidate}"
            break
        fi
    done
fi

if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "Error: executable not found or not executable: ${EXECUTABLE}" >&2
    echo "Build it with: make -C \"${ROOT_DIR}/exec\" -f GNUmakefile_TS -j4" >&2
    exit 1
fi
if [[ ! -f "${BASE_INPUT}" ]]; then
    echo "Error: base inputs file not found: ${BASE_INPUT}" >&2
    exit 1
fi
if [[ ! -f "${PLOTTER}" ]]; then
    echo "Error: plotting script not found: ${PLOTTER}" >&2
    exit 1
fi
if ! command -v "${MPIEXEC}" >/dev/null 2>&1; then
    echo "Error: MPI launcher not found: ${MPIEXEC}" >&2
    exit 1
fi
if [[ -z "${PYTHON_BIN}" ]] || ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    echo "Error: no Python interpreter with NumPy and Matplotlib was found." >&2
    echo "Create the project plotting environment with:" >&2
    echo "  python3 -m venv \"${ROOT_DIR}/.venv\"" >&2
    echo "  \"${ROOT_DIR}/.venv/bin/python\" -m pip install -r \"${SCRIPT_DIR}/requirements-plotting.txt\"" >&2
    echo "Or set PYTHON_BIN to an interpreter that provides those packages." >&2
    exit 1
fi
if ! "${PYTHON_BIN}" -c 'import matplotlib, numpy' >/dev/null 2>&1; then
    echo "Error: ${PYTHON_BIN} does not provide NumPy and Matplotlib." >&2
    echo "Install them with:" >&2
    echo "  \"${PYTHON_BIN}\" -m pip install -r \"${SCRIPT_DIR}/requirements-plotting.txt\"" >&2
    exit 1
fi
mkdir -p "${RESULTS_DIR}"

echo "Selected tests: ${TEST_CASES[*]}"
echo "Mesh: domainSize=${DOMAIN_SIZE}, thickness=${THICKNESS}"
echo "Boxes: boxSize_nonrad=${BOX_SIZE_NONRAD}, boxSize_rad=${BOX_SIZE_RAD}"
echo "MPI processes: ${NPROCS}"

stream_simulation_progress() {
    local case_number="$1"
    local case_name="$2"
    local last_iteration="-1"
    local line
    local iteration
    local percent

    while IFS= read -r line; do
        printf '%s\n' "${line}"
        if [[ "${line}" =~ iter[[:space:]]*=[[:space:]]*([0-9]+) ]]; then
            iteration="${BASH_REMATCH[1]}"
            if [[ "${iteration}" != "${last_iteration}" ]]; then
                percent=$((iteration * 100 / MAX_ITER))
                if (( percent > 100 )); then
                    percent=100
                fi
                printf '[case %d/%d: %s] simulation %d/%d (%d%%)\n' \
                    "${case_number}" "${TOTAL_CASES}" "${case_name}" \
                    "${iteration}" "${MAX_ITER}" "${percent}"
                last_iteration="${iteration}"
            fi
        fi
    done
}

make_case_input() {
    local case_id="$1"
    local case_name="$2"
    local output_file="$3"
    local max_time="$4"

    awk \
        -v case_id="${case_id}" \
        -v data_prefix="${case_name}" \
        -v checkpoint_prefix="${case_name}_checkpoint" \
        -v domain_size="${DOMAIN_SIZE}" \
        -v thickness="${THICKNESS}" \
        -v box_size_nonrad="${BOX_SIZE_NONRAD}" \
        -v box_size_rad="${BOX_SIZE_RAD}" \
        -v max_iter="${MAX_ITER}" \
        -v max_time="${max_time}" \
        -v slice_cadence="${SLICE_CADENCE}" \
        '
        BEGIN { found_slices = 0 }
        $1 == "-init_condition_type"  { $2 = case_id }
        $1 == "-gamma" && (case_id == 8 || case_id == 9) { $2 = "1.4" }
        $1 == "-max_time"             { $2 = max_time }
        $1 == "-domainSize"           { $2 = domain_size }
        $1 == "-thickness"            { $2 = thickness }
        $1 == "-boxSize_nonrad"       { $2 = box_size_nonrad }
        $1 == "-boxSize_rad"          { $2 = box_size_rad }
        $1 == "-max_iter"             { $2 = max_iter }
        $1 == "-slice_cadence"        { $2 = slice_cadence }
        $1 == "-slice_time_cadence"   { $2 = "1.0e30" }
        $1 == "-slices"               { $2 = "Z"; found_slices = 1 }
        $1 == "-write_cadence"        { $2 = max_iter }
        $1 == "-write_time_cadence"   { $2 = "1.0e30" }
        $1 == "-checkpoint_cadence"   { $2 = max_iter }
        $1 == "-data_file_prefix"     { $2 = data_prefix }
        $1 == "-checkpoint_file_prefix" { $2 = checkpoint_prefix }
        { print }
        END {
            if (!found_slices) print "-slices Z"
        }
        ' "${BASE_INPUT}" > "${output_file}"
}

clean_case_directory() {
    local case_dir="$1"
    local case_name="$2"
    local expected_case_dir="${RESULTS_DIR}/${case_name}"

    # case_name comes from the fixed mapping above. Keep this guard close to
    # the recursive deletion so a future path change cannot broaden its scope.
    if [[ -z "${case_dir}" || "${case_dir}" == "/" \
        || "${case_dir}" == "${RESULTS_DIR}" \
        || "${case_dir}" != "${expected_case_dir}" ]]; then
        echo "Error: refusing to clean unsafe case directory: '${case_dir}'." >&2
        exit 1
    fi

    if [[ -e "${case_dir}" || -L "${case_dir}" ]]; then
        echo "Cleaning previous results: ${case_dir}"
        rm -rf -- "${case_dir}"
    fi
    mkdir -p "${case_dir}"
}

REGULAR_MAX_ITER="${MAX_ITER}"
for index in "${!CASE_IDS[@]}"; do
    case_number=$((index + 1))
    case_id="${CASE_IDS[${index}]}"
    MAX_ITER="${REGULAR_MAX_ITER}"
    if [[ "${case_id}" == 8 || "${case_id}" == 9 ]]; then
        MAX_ITER="${BLAST_MAX_ITER}"
    fi
    case_name="${CASE_NAMES[${index}]}"
    case_max_time="$(max_time_for_id "${case_id}")"
    case_dir="${RESULTS_DIR}/${case_name}"
    case_input="${case_dir}/inputs"
    plot_dir="${case_dir}/plots_z0_raw"

    clean_case_directory "${case_dir}" "${case_name}"
    make_case_input "${case_id}" "${case_name}" "${case_input}" "${case_max_time}"

    echo
    echo "=== [case ${case_number}/${TOTAL_CASES}] Running ${case_name} (init_condition_type=${case_id}) ==="
    echo "Results: ${case_dir}"
    echo "Limits: max_iter=${MAX_ITER}, max_time=${case_max_time}"
    case_start="${SECONDS}"
    (
        cd "${case_dir}"
        "${MPIEXEC}" -np "${NPROCS}" "${EXECUTABLE}" "${case_input}" \
            2>&1 \
            | stream_simulation_progress "${case_number}" "${case_name}" \
            | tee run.log
    )

    echo "=== [case ${case_number}/${TOTAL_CASES}] Plotting raw z=0 slices for ${case_name} ==="
    "${PYTHON_BIN}" "${PLOTTER}" "${case_dir}" \
        --output-dir "${plot_dir}" \
        --workers "${PLOT_WORKERS}"
    case_elapsed=$((SECONDS - case_start))
    echo "=== [case ${case_number}/${TOTAL_CASES}] Completed ${case_name} in ${case_elapsed}s ==="
done

echo
echo "All ${TOTAL_CASES} selected test problem(s) completed."
echo "Results and raw z=0 plots are under: ${RESULTS_DIR}"
