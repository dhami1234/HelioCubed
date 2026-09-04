#!/usr/bin/env bash

# Run the four idealized HelioCubed test problems and plot their z=0 slices.
#
# Usage:
#   ./scripts/run_test_problems.sh [results_directory]
#
# Defaults:
#   - Runs init_condition_type 0, 1, 5, and 6 sequentially.
#   - Uses 6 MPI processes, 100 iterations, and writes z=0 slices every
#     10 iterations (including the initial and final states).
#   - Writes each case into its own directory under ./Test_results.
#   - Plots density, velocity, pressure, magnetic field, and derived
#     temperature in raw solver/CGS units under each case's plots_z0_raw/.
#
# Resolution:
#   The script does not set the mesh resolution. Every case inherits these
#   values from exec/inputs_Test_Problems:
#
#     domainSize       cells in each angular direction on every cubed-sphere face
#     thickness        cells in the radial direction
#     boxSize_nonrad   angular cells per computational/MPI box
#     boxSize_rad      radial cells per computational/MPI box
#
#   The approximate cell count is:
#
#     6 * domainSize * domainSize * thickness
#
#   boxSize_nonrad and boxSize_rad change domain decomposition, not physical
#   resolution. Their values should divide domainSize and thickness cleanly.
#   NPROCS changes parallelism only; it does not change the mesh resolution.
#
# Useful overrides:
#   NPROCS=8 PLOT_WORKERS=4 ./scripts/run_test_problems.sh
#   PYTHON_BIN=/path/to/python ./scripts/run_test_problems.sh
#   HELIOCUBED_EXE=/path/to/exe ./scripts/run_test_problems.sh /path/to/results
#   HELIOCUBED_TEST_INPUT=/path/to/inputs ./scripts/run_test_problems.sh
#
# The selected Python must provide NumPy and Matplotlib. Existing result
# directories are reused: generated inputs, logs, and matching output files
# are overwritten, but unrelated or stale files are not deleted automatically.
# The script stops immediately if a simulation or plotting step fails.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
EXECUTABLE="${HELIOCUBED_EXE:-${ROOT_DIR}/exec/cubedSphereTest.exe}"
BASE_INPUT="${HELIOCUBED_TEST_INPUT:-${ROOT_DIR}/exec/inputs_Test_Problems}"
RESULTS_DIR="${1:-${ROOT_DIR}/Test_results}"
MPIEXEC="${MPIEXEC:-mpirun}"
NPROCS="${NPROCS:-6}"
PLOT_WORKERS="${PLOT_WORKERS:-8}"
PLOTTER="${SCRIPT_DIR}/Plot_Test_Slices.py"

if [[ -n "${PYTHON_BIN:-}" ]]; then
    PYTHON_BIN="${PYTHON_BIN}"
else
    PYTHON_BIN=""
    for python_candidate in python3 /usr/local/bin/python3.10; do
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
    echo "Set PYTHON_BIN to the interpreter used for scripts/Plot_Slices.py." >&2
    exit 1
fi
if ! "${PYTHON_BIN}" -c 'import matplotlib, numpy' >/dev/null 2>&1; then
    echo "Error: ${PYTHON_BIN} does not provide NumPy and Matplotlib." >&2
    echo "Set PYTHON_BIN to the interpreter used for scripts/Plot_Slices.py." >&2
    exit 1
fi
if ! [[ "${NPROCS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: NPROCS must be a positive integer (got '${NPROCS}')." >&2
    exit 1
fi
if ! [[ "${PLOT_WORKERS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: PLOT_WORKERS must be a positive integer (got '${PLOT_WORKERS}')." >&2
    exit 1
fi

mkdir -p "${RESULTS_DIR}"

CASE_IDS=(0 1 5 6)
CASE_NAMES=(
    smooth_spherical_wave
    non_radial_spherical_wave
    symmetric_strong_shock
    strong_shock_constant_B
)
TOTAL_CASES="${#CASE_IDS[@]}"

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
                percent=$((iteration * 100 / 100))
                printf '[case %d/%d: %s] simulation %d/100 (%d%%)\n' \
                    "${case_number}" "${TOTAL_CASES}" "${case_name}" \
                    "${iteration}" "${percent}"
                last_iteration="${iteration}"
            fi
        fi
    done
}

make_case_input() {
    local case_id="$1"
    local case_name="$2"
    local output_file="$3"

    awk \
        -v case_id="${case_id}" \
        -v data_prefix="${case_name}" \
        -v checkpoint_prefix="${case_name}_checkpoint" \
        '
        BEGIN { found_slices = 0 }
        $1 == "-init_condition_type"  { $2 = case_id }
        $1 == "-max_iter"             { $2 = 100 }
        $1 == "-slice_cadence"        { $2 = 10 }
        $1 == "-slice_time_cadence"   { $2 = "1.0e30" }
        $1 == "-slices"               { $2 = "Z"; found_slices = 1 }
        $1 == "-write_cadence"        { $2 = 100 }
        $1 == "-write_time_cadence"   { $2 = "1.0e30" }
        $1 == "-checkpoint_cadence"   { $2 = 100 }
        $1 == "-data_file_prefix"     { $2 = data_prefix }
        $1 == "-checkpoint_file_prefix" { $2 = checkpoint_prefix }
        { print }
        END {
            if (!found_slices) print "-slices Z"
        }
        ' "${BASE_INPUT}" > "${output_file}"
}

for index in "${!CASE_IDS[@]}"; do
    case_number=$((index + 1))
    case_id="${CASE_IDS[${index}]}"
    case_name="${CASE_NAMES[${index}]}"
    case_dir="${RESULTS_DIR}/${case_name}"
    case_input="${case_dir}/inputs"
    plot_dir="${case_dir}/plots_z0_raw"

    mkdir -p "${case_dir}"
    make_case_input "${case_id}" "${case_name}" "${case_input}"

    echo
    echo "=== [case ${case_number}/${TOTAL_CASES}] Running ${case_name} (init_condition_type=${case_id}) ==="
    echo "Results: ${case_dir}"
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
echo "All four test problems completed."
echo "Results and raw z=0 plots are under: ${RESULTS_DIR}"
