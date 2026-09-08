#!/usr/bin/env python3
"""Small solver integration tests; build the executable before running this file.

Uses isolated temporary output directories and a single MPI process. Override
HELIOCUBED_EXE to test another build. No Python packages beyond the standard
library are required.
"""

import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(os.environ.get("HELIOCUBED_EXE", ROOT / "exec/cubedSphereTest.exe")).resolve()
STEP = 2.0 ** -13  # Exactly representable even if the input parser uses floats.


class ConvergenceTimeLimitTest(unittest.TestCase):
    def run_case(self, **overrides):
        temporary = tempfile.TemporaryDirectory(prefix="hc-time-limit-")
        self.addCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        settings = dict(
            convTestType=1, convergence_level=0, convergence_dt=STEP,
            init_condition_type=2, domainSize=8, thickness=8,
            boxSize_nonrad=8, boxSize_rad=8, max_iter=4,
            max_time=1.5 * STEP, temporal_order=1,
            write_cadence=1000, slice_cadence=1000, checkpoint_cadence=1000,
            write_time_cadence=1.e30, slice_time_cadence=1.e30,
            P_floor_cadence=1000, slices="Z",
        )
        settings.update(overrides)
        lines = []
        for line in (ROOT / "exec/inputs_convergence").read_text().splitlines():
            fields = line.split()
            if fields and fields[0].startswith("-"):
                key = fields[0][1:]
                if key in settings:
                    line = f"-{key} {settings.pop(key)}"
            lines.append(line)
        lines.extend(f"-{key} {value}" for key, value in settings.items())
        (directory / "inputs").write_text("\n".join(lines) + "\n")
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        completed = subprocess.run(
            [str(EXE), "inputs"], cwd=directory, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            timeout=180,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout[-6000:])
        log = completed.stdout
        steps = [
            (int(i), float(dt), float(t))
            for i, dt, t in re.findall(
                r"iter = (\d+) dt = ([\deE+.-]+)\(s\) time = ([\deE+.-]+)\(s\)", log
            )
        ]
        return directory, log, steps

    def assert_final_files(self, directory, iteration):
        for name in (f"data.{iteration}.hdf5", f"data.z.{iteration}.dat",
                     f"checkpoint_{iteration}.hdf5", "U_conv_test_0.hdf5"):
            self.assertTrue((directory / name).is_file(), name)

    def test_fractional_last_step_and_final_output(self):
        directory, log, steps = self.run_case()
        self.assertEqual([s[0] for s in steps], [1, 2])
        self.assertAlmostEqual(steps[-1][1], STEP / 2, delta=STEP * 1.e-6)
        self.assertIn("Maximum simulation time reached", log)
        self.assert_final_files(directory, 2)

    def test_exact_boundary_does_not_take_zero_step(self):
        directory, log, steps = self.run_case(max_time=2 * STEP)
        self.assertEqual([s[0] for s in steps], [1, 2])
        self.assertTrue(all(s[1] > 0 for s in steps))
        self.assertIn("Maximum simulation time reached", log)
        self.assert_final_files(directory, 2)

    def test_iteration_limit_still_applies(self):
        directory, log, steps = self.run_case(max_iter=1, max_time=10 * STEP)
        self.assertEqual(len(steps), 1)
        self.assertAlmostEqual(steps[0][1], STEP, delta=STEP * 1.e-6)
        self.assertNotIn("Maximum simulation time reached", log)
        self.assert_final_files(directory, 1)

    def test_restart_at_time_limit_does_not_advance(self):
        directory, _, _ = self.run_case()
        _, _, steps = self.run_case(restart_file=directory / "checkpoint_2.hdf5")
        self.assertEqual(steps, [])

    def test_nominal_timestep_preserved_between_levels(self):
        for mode in (1, 2):
            with self.subTest(mode=mode):
                directory, log, steps = self.run_case(
                    convTestType=mode, convergence_level=-1, convergence_dt=-1,
                    max_time=STEP / 2,
                )
                timesteps = [float(dt) for dt in re.findall(
                    r"Convergence timestep = ([\deE+.-]+)", log
                )]
                self.assertEqual(len(timesteps), 3)
                self.assertGreater(timesteps[0], STEP / 2)
                for level, dt in enumerate(timesteps):
                    expected = timesteps[0] / (2 ** level if mode == 2 else 1)
                    self.assertTrue(math.isclose(dt, expected, rel_tol=1.e-12))
                self.assertEqual(len(steps), 3)
                self.assertEqual(log.count("Maximum simulation time reached"), 3)
                self.assertTrue(all(math.isclose(s[2], STEP / 2, rel_tol=5.e-4)
                                    for s in steps))
                for level in range(3):
                    self.assertTrue((directory / f"U_conv_test_{level}.hdf5").is_file())


if __name__ == "__main__":
    unittest.main()
