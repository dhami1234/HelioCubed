#!/bin/bash

#SBATCH --nodes=7
#SBATCH --ntasks=864
#SBATCH --time=24:00:00
#SBATCH --job-name=heliocubed-test
#SBATCH --output=myjob.o%j
#SBATCH --error=myjob.e%j
#SBATCH --partition=wholenode

# Submit selected idealized HelioCubed test problems as separate Anvil jobs.
#
# Run this script from the directory that should contain the case directories:
#
#   cd /path/to/test/results
#   /anvil/projects/x-mca07s033/Talwinder/HelioCubed/scripts/run_test_problems_anvil.sh
#
# Each selected case gets its own directory, generated input file, and Slurm
# job. The same script is submitted to Slurm in worker mode; do not invoke
# --run-case manually.

set -Eeuo pipefail

run_anvil_case() {
    local executable="$1"
    local input_file="$2"

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
    command -v srun

    export MV2_HOMOGENEOUS_CLUSTER=1
    export OMP_NUM_THREADS=1
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close

    # Shared libraries used by the Anvil build. Environment overrides make it
    # possible to update installation paths without editing this function.
    local hdf5_lib="${HELIOCUBED_ANVIL_HDF5_LIB:-/apps/spack/anvil/apps/hdf5/1.10.7-intel-19.0.5-gddu73a/lib}"
    local gsl_lib="${HELIOCUBED_ANVIL_GSL_LIB:-/home/x-dhami/gsl/lib}"
    local openblas_lib="${HELIOCUBED_ANVIL_OPENBLAS_LIB:-/apps/spack/anvil/apps/openblas/0.3.17-gcc-11.2.0-2qrsari/lib}"

    if [[ ! -r "${openblas_lib}/libopenblas.so.0" ]]; then
        echo "ERROR: Cannot find ${openblas_lib}/libopenblas.so.0" >&2
        echo "Check the available installation with: module spider openblas/0.3.17" >&2
        exit 1
    fi
    export LD_LIBRARY_PATH="${openblas_lib}:${gsl_lib}:${hdf5_lib}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

    if [[ ! -x "${executable}" ]]; then
        echo "ERROR: executable not found or not executable: ${executable}" >&2
        exit 1
    fi
    if [[ ! -f "${input_file}" ]]; then
        echo "ERROR: input file not found: ${input_file}" >&2
        exit 1
    fi

    cd "${SLURM_SUBMIT_DIR:-$(dirname "${input_file}")}"
    echo "Starting ${SLURM_JOB_NAME:-HelioCubed test} on $(date)"
    echo "Working directory: $(pwd)"
    echo "MPI tasks: ${SLURM_NTASKS:-unknown}"

    srun --mpi=pmi2 -n "${SLURM_NTASKS}" "${executable}" "${input_file}"

    echo "Finished on $(date)"
}

if [[ "${1:-}" == "--run-case" ]]; then
    if [[ "$#" -ne 3 ]]; then
        echo "ERROR: internal worker invocation requires an executable and input file." >&2
        exit 2
    fi
    run_anvil_case "$2" "$3"
    exit 0
fi

# =============================================================================
# USER CONFIGURATION -- edit these values before submitting
# =============================================================================

# Anvil directory containing cubedSphereTest.exe and inputs_Test_Problems.
EXEC_DIR="${HELIOCUBED_EXEC_DIR:-/anvil/projects/x-mca07s033/Talwinder/HelioCubed/exec}"

# Tests to run. Use one or more IDs, for example TEST_CASES=(4) or (2 3 8).
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

# Mesh resolution and computational/MPI box dimensions. The box sizes must
# divide the corresponding mesh dimensions exactly.
DOMAIN_SIZE="${HELIOCUBED_TEST_DOMAIN_SIZE:-60}"
THICKNESS="${HELIOCUBED_TEST_THICKNESS:-120}"
BOX_SIZE_NONRAD="${HELIOCUBED_TEST_BOX_SIZE_NONRAD:-60}"
BOX_SIZE_RAD="${HELIOCUBED_TEST_BOX_SIZE_RAD:-40}"

# Simulation limits and output cadence.
MAX_ITER="${HELIOCUBED_TEST_MAX_ITER:-100}"
BLAST_MAX_ITER="${HELIOCUBED_TEST_BLAST_MAX_ITER:-10000}"
SLICE_CADENCE="${HELIOCUBED_TEST_SLICE_CADENCE:-50}"

# Maximum simulation time for each test problem.
CASE_2_MAX_TIME="${HELIOCUBED_TEST_CASE_2_MAX_TIME:-3600000.0}" # radial pulse
CASE_3_MAX_TIME="${HELIOCUBED_TEST_CASE_3_MAX_TIME:-3600000.0}" # radial pulse, constant B
CASE_4_MAX_TIME="${HELIOCUBED_TEST_CASE_4_MAX_TIME:-3600000.0}" # non-radial pulse
CASE_5_MAX_TIME="${HELIOCUBED_TEST_CASE_5_MAX_TIME:-3600000.0}" # non-radial pulse, constant B
CASE_6_MAX_TIME="${HELIOCUBED_TEST_CASE_6_MAX_TIME:-3600000.0}" # strong spherical shock
CASE_7_MAX_TIME="${HELIOCUBED_TEST_CASE_7_MAX_TIME:-3600000.0}" # strong shock, constant B
CASE_8_MAX_TIME="${HELIOCUBED_TEST_CASE_8_MAX_TIME:-${HELIOCUBED_TEST_BLAST_MAX_TIME:-5001.23}}" # hydro blast
CASE_9_MAX_TIME="${HELIOCUBED_TEST_CASE_9_MAX_TIME:-${HELIOCUBED_TEST_BLAST_MAX_TIME:-5001.23}}" # MHD blast

# Slurm resources. Command-line sbatch options generated below override the
# matching #SBATCH defaults at the top of this file.
NODES="${HELIOCUBED_ANVIL_NODES:-7}"
# NPROCS is passed to Slurm as --ntasks and is the number of MPI processes.
NPROCS="${HELIOCUBED_TEST_NPROCS:-${HELIOCUBED_ANVIL_NTASKS:-864}}"
WALLTIME="${HELIOCUBED_ANVIL_WALLTIME:-24:00:00}"
PARTITION="${HELIOCUBED_ANVIL_PARTITION:-wholenode}"

# Like the laptop runner, remove a selected case's existing directory before
# recreating it. Set to 0 to stop instead when a selected directory exists.
CLEAN_EXISTING_RESULTS="${HELIOCUBED_TEST_CLEAN_EXISTING_RESULTS:-1}"

# =============================================================================
# END USER CONFIGURATION
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
EXECUTABLE="${HELIOCUBED_EXE:-${EXEC_DIR}/cubedSphereTest.exe}"
BASE_INPUT="${HELIOCUBED_TEST_INPUT:-${EXEC_DIR}/inputs_Test_Problems}"
RESULTS_DIR="$(pwd -P)"

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
        echo "ERROR: ${variable_name} must be a positive integer (got '${value}')." >&2
        exit 1
    fi
}

require_positive_number() {
    local variable_name="$1"
    local value="$2"

    if ! awk -v value="${value}" 'BEGIN {
        exit !(value ~ /^[0-9]+([.][0-9]*)?([eE][+-]?[0-9]+)?$/ && value+0 > 0)
    }'; then
        echo "ERROR: ${variable_name} must be a positive number (got '${value}')." >&2
        exit 1
    fi
}

require_walltime() {
    if ! [[ "$1" =~ ^[0-9]+:[0-5][0-9]:[0-5][0-9]$ ]]; then
        echo "ERROR: WALLTIME must use hh:mm:ss format (got '$1')." >&2
        exit 1
    fi
}

make_case_input() {
    local case_id="$1"
    local case_name="$2"
    local output_file="$3"
    local max_iter="$4"
    local max_time="$5"

    awk \
        -v case_id="${case_id}" \
        -v data_prefix="${case_name}" \
        -v checkpoint_prefix="${case_name}_checkpoint" \
        -v domain_size="${DOMAIN_SIZE}" \
        -v thickness="${THICKNESS}" \
        -v box_size_nonrad="${BOX_SIZE_NONRAD}" \
        -v box_size_rad="${BOX_SIZE_RAD}" \
        -v max_iter="${max_iter}" \
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

prepare_case_directory() {
    local case_dir="$1"
    local case_name="$2"
    local expected_case_dir="${RESULTS_DIR}/${case_name}"

    if [[ -z "${case_dir}" || "${case_dir}" == "/" \
        || "${case_dir}" == "${RESULTS_DIR}" \
        || "${case_dir}" != "${expected_case_dir}" ]]; then
        echo "ERROR: refusing to clean unsafe case directory: '${case_dir}'." >&2
        exit 1
    fi

    if [[ -e "${case_dir}" || -L "${case_dir}" ]]; then
        if [[ "${CLEAN_EXISTING_RESULTS}" != 1 ]]; then
            echo "ERROR: selected case directory already exists: ${case_dir}" >&2
            echo "Remove it or set CLEAN_EXISTING_RESULTS=1." >&2
            exit 1
        fi
        echo "Cleaning previous results: ${case_dir}"
        rm -rf -- "${case_dir}"
    fi
    mkdir -p "${case_dir}"
}

if [[ "${#TEST_CASES[@]}" -eq 0 ]]; then
    echo "ERROR: TEST_CASES must contain at least one test ID." >&2
    exit 1
fi
if [[ "${CLEAN_EXISTING_RESULTS}" != 0 && "${CLEAN_EXISTING_RESULTS}" != 1 ]]; then
    echo "ERROR: CLEAN_EXISTING_RESULTS must be 0 or 1." >&2
    exit 1
fi

require_positive_integer "DOMAIN_SIZE" "${DOMAIN_SIZE}"
require_positive_integer "THICKNESS" "${THICKNESS}"
require_positive_integer "BOX_SIZE_NONRAD" "${BOX_SIZE_NONRAD}"
require_positive_integer "BOX_SIZE_RAD" "${BOX_SIZE_RAD}"
require_positive_integer "MAX_ITER" "${MAX_ITER}"
require_positive_integer "BLAST_MAX_ITER" "${BLAST_MAX_ITER}"
require_positive_integer "SLICE_CADENCE" "${SLICE_CADENCE}"
require_positive_integer "NODES" "${NODES}"
require_positive_integer "NPROCS" "${NPROCS}"
require_walltime "${WALLTIME}"

require_positive_number "CASE_2_MAX_TIME" "${CASE_2_MAX_TIME}"
require_positive_number "CASE_3_MAX_TIME" "${CASE_3_MAX_TIME}"
require_positive_number "CASE_4_MAX_TIME" "${CASE_4_MAX_TIME}"
require_positive_number "CASE_5_MAX_TIME" "${CASE_5_MAX_TIME}"
require_positive_number "CASE_6_MAX_TIME" "${CASE_6_MAX_TIME}"
require_positive_number "CASE_7_MAX_TIME" "${CASE_7_MAX_TIME}"
require_positive_number "CASE_8_MAX_TIME" "${CASE_8_MAX_TIME}"
require_positive_number "CASE_9_MAX_TIME" "${CASE_9_MAX_TIME}"

if (( DOMAIN_SIZE % BOX_SIZE_NONRAD != 0 )); then
    echo "ERROR: BOX_SIZE_NONRAD (${BOX_SIZE_NONRAD}) must divide DOMAIN_SIZE (${DOMAIN_SIZE})." >&2
    exit 1
fi
if (( THICKNESS % BOX_SIZE_RAD != 0 )); then
    echo "ERROR: BOX_SIZE_RAD (${BOX_SIZE_RAD}) must divide THICKNESS (${THICKNESS})." >&2
    exit 1
fi
if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "ERROR: executable not found or not executable: ${EXECUTABLE}" >&2
    exit 1
fi
if [[ ! -f "${BASE_INPUT}" ]]; then
    echo "ERROR: base input file not found: ${BASE_INPUT}" >&2
    exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is not available. Run this launcher on an Anvil login node." >&2
    exit 1
fi

CASE_IDS=("${TEST_CASES[@]}")
CASE_NAMES=()
seen_case_ids=" "
for case_id in "${CASE_IDS[@]}"; do
    if ! case_name="$(case_name_for_id "${case_id}")"; then
        echo "ERROR: unsupported test ID '${case_id}'; choose from 2 through 9." >&2
        exit 1
    fi
    if [[ "${seen_case_ids}" == *" ${case_id} "* ]]; then
        echo "ERROR: test ID '${case_id}' appears more than once in TEST_CASES." >&2
        exit 1
    fi
    seen_case_ids+="${case_id} "
    CASE_NAMES+=("${case_name}")
done

echo "Submitting tests from: ${RESULTS_DIR}"
echo "Selected tests: ${TEST_CASES[*]}"
echo "Mesh: domainSize=${DOMAIN_SIZE}, thickness=${THICKNESS}"
echo "Boxes: boxSize_nonrad=${BOX_SIZE_NONRAD}, boxSize_rad=${BOX_SIZE_RAD}"
echo "Slurm: nodes=${NODES}, MPI processes=${NPROCS}, time=${WALLTIME}, partition=${PARTITION}"

submitted_count=0
for index in "${!CASE_IDS[@]}"; do
    case_id="${CASE_IDS[${index}]}"
    case_name="${CASE_NAMES[${index}]}"
    case_max_iter="${MAX_ITER}"
    if [[ "${case_id}" == 8 || "${case_id}" == 9 ]]; then
        case_max_iter="${BLAST_MAX_ITER}"
    fi
    case_max_time="$(max_time_for_id "${case_id}")"
    case_dir="${RESULTS_DIR}/${case_name}"
    case_input="${case_dir}/inputs"
    job_name="hc_t${case_id}"

    prepare_case_directory "${case_dir}" "${case_name}"
    make_case_input \
        "${case_id}" "${case_name}" "${case_input}" \
        "${case_max_iter}" "${case_max_time}"

    submission="$(
        cd "${case_dir}"
        sbatch \
            --nodes="${NODES}" \
            --ntasks="${NPROCS}" \
            --time="${WALLTIME}" \
            --partition="${PARTITION}" \
            --job-name="${job_name}" \
            --output="${case_dir}/myjob.o%j" \
            --error="${case_dir}/myjob.e%j" \
            --chdir="${case_dir}" \
            "${SCRIPT_PATH}" --run-case "${EXECUTABLE}" "${case_input}"
    )"
    submitted_count=$((submitted_count + 1))

    echo "[${submitted_count}/${#CASE_IDS[@]}] ${case_name}: ${submission}"
    echo "    max_iter=${case_max_iter}, max_time=${case_max_time}"
    echo "    directory=${case_dir}"
done

echo "Submitted ${submitted_count} test problem job(s)."
echo "Use 'squeue -u ${USER:-$(id -un)}' to monitor them."
