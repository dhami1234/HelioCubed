#!/usr/bin/env python3
"""Plot HelioCubed z=0 test-problem slices without output scaling.

This is the test-problem counterpart to Plot_Slices.py.  It reuses that
script's Tecplot parser and cubed-sphere plotting geometry, while deliberately
disabling unit conversion, radial scaling, and corotation so the color bars
show the values written by HelioCubed in its native CGS units.
"""

from __future__ import annotations

import argparse
import os
from multiprocessing import Pool, cpu_count
from pathlib import Path
import re
import sys

# The default ~/.matplotlib is not writable on some compute nodes.  Use a
# per-user temporary cache unless the caller has selected another location.
matplotlib_cache = Path(os.environ.get("TMPDIR", "/tmp")) / (
    f"heliocubed-matplotlib-{os.getuid()}"
)
matplotlib_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))

import matplotlib

matplotlib.use("Agg")

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import Plot_Slices as base_plotter


VARIABLES = ["density", "Vr", "Vt", "Vp", "P", "Br", "Bt", "Bp", "T"]

RAW_LABELS = {
    "density": "density [g/cm^3]",
    "Vr": "Vr [cm/s]",
    "Vt": "Vt [cm/s]",
    "Vp": "Vp [cm/s]",
    "P": "P [dyn/cm^2]",
    "Br": "Br [G]",
    "Bt": "Bt [G]",
    "Bp": "Bp [G]",
    "T": "T [K]",
}

RAW_COLORMAPS = {
    "density": {"cmap": "terrain_r"},
    "Vr": {"cmap": "rainbow"},
    "Vt": {"cmap": "bwr"},
    "Vp": {"cmap": "bwr"},
    "P": {"cmap": "plasma"},
    "Br": {"cmap": "bwr"},
    "Bt": {"cmap": "bwr"},
    "Bp": {"cmap": "bwr"},
    "T": {"cmap": "hot"},
}


def configure_raw_output() -> None:
    """Disable every value scaling operation in Plot_Slices.py."""
    base_plotter.COROTATE = False
    base_plotter.UNIT_SCALE = {name: 1.0 for name in VARIABLES}
    base_plotter.RADIAL_SCALE = {name: 0.0 for name in VARIABLES}
    base_plotter.LABELS = RAW_LABELS
    # Omitting vmin/vmax lets each raw data set choose an appropriate range.
    base_plotter.COLOR_SETTINGS = RAW_COLORMAPS


def render_slice(task: tuple[str, str, tuple[str, ...]]) -> list[str]:
    """Render every requested variable from one z=0 slice file."""
    filepath, output_dir, variables = task
    configure_raw_output()
    geometry = base_plotter.prepare_geometry(
        filename=filepath,
        reshape_order="C",
        corotate=False,
    )
    return [
        base_plotter.plot_with_geometry(
            geometry,
            varname=variable,
            cell_order="C",
            out_path=output_dir,
        )
        for variable in variables
    ]


def iteration_number(path: Path) -> int:
    match = re.search(r"\.z\.(\d+)\.dat$", path.name)
    if match is None:
        raise ValueError(f"Not a z=0 slice filename: {path}")
    return int(match.group(1))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot all raw variables from HelioCubed z=0 slices."
    )
    parser.add_argument("slice_directory", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="PNG directory (default: <slice_directory>/plots_z0_raw)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=min(8, max(1, (cpu_count() or 2) - 1)),
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        choices=VARIABLES,
        default=VARIABLES,
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    slice_directory = args.slice_directory.resolve()
    output_directory = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else slice_directory / "plots_z0_raw"
    )
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    if not slice_directory.is_dir():
        raise FileNotFoundError(f"Slice directory does not exist: {slice_directory}")

    slice_files = sorted(
        slice_directory.glob("*.z.*.dat"), key=iteration_number
    )
    if not slice_files:
        raise FileNotFoundError(
            f"No z=0 slices matching '*.z.*.dat' in {slice_directory}"
        )

    output_directory.mkdir(parents=True, exist_ok=True)
    variables = tuple(args.variables)
    tasks = [
        (str(slice_file), str(output_directory), variables)
        for slice_file in slice_files
    ]

    total = len(tasks) * len(variables)
    print(
        f"[INFO] Plotting {len(tasks)} z=0 slices x {len(variables)} variables "
        f"in raw solver units ({total} PNGs)."
    )

    completed = 0
    if args.workers == 1:
        results = map(render_slice, tasks)
        for outputs in results:
            completed += len(outputs)
            print(f"[{completed}/{total}] {Path(outputs[0]).stem.rsplit('_', 1)[-1]}")
    else:
        with Pool(processes=args.workers) as pool:
            for outputs in pool.imap_unordered(render_slice, tasks, chunksize=1):
                completed += len(outputs)
                print(f"[{completed}/{total}] {Path(outputs[0]).stem.rsplit('_', 1)[-1]}")

    print(f"[INFO] Raw z=0 plots written to {output_directory}")


if __name__ == "__main__":
    main()
