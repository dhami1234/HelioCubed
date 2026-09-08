#!/usr/bin/env python3
"""Solver integration checks for cases 8/9; build the executable first.

Runs in temporary directories with one MPI rank and reads the raw z=0 slices.
No Python packages beyond the standard library are required. The expected
physical values assume the default shell geometry and MHD_Constants.H values.
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


def read_slice(path):
    text = path.read_text()
    names = re.findall(r'"([^"]+)"', re.search(r'^VARIABLES=(.*)$', text, re.M)[1])[3:]
    cells = {name: [] for name in names}
    for zone in re.split(r'^ZONE ', text, flags=re.M)[1:]:
        header, data = zone.split('\n', 1)
        ni, nj = (int(re.search(r'\b' + axis + r'=(\d+)', header)[1]) for axis in ('I', 'J'))
        values = [float(v) for v in data.split()]
        n = (ni - 1) * (nj - 1)
        assert len(values) == 3 * ni * nj + len(names) * n
        offset = 3 * ni * nj
        for name in names:
            cells[name].extend(values[offset:offset+n])
            offset += n
    return cells


class BlastProblemsTest(unittest.TestCase):
    def run_case(self, case):
        temporary = tempfile.TemporaryDirectory(prefix=f'hc-blast-{case}-')
        self.addCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        source = (ROOT / 'exec/inputs_Blast').read_text()
        settings = dict(init_condition_type=case, domainSize=16, thickness=16,
                        boxSize_nonrad=16, boxSize_rad=16, max_time=10,
                        max_iter=100, slice_cadence=1)
        for key, value in settings.items():
            source = re.sub(r'^-' + key + r'\s+.*$', f'-{key} {value}', source, flags=re.M)
        (directory / 'inputs').write_text(source)
        completed = subprocess.run(
            [str(EXE), 'inputs'], cwd=directory,
            env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=180)
        self.assertEqual(completed.returncode, 0, completed.stdout[-6000:])
        self.assertIn('Maximum simulation time reached', completed.stdout)
        files = sorted(directory.glob('blast_data.z.*.dat'), key=lambda p: int(p.name.split('.')[-2]))
        self.assertGreaterEqual(len(files), 2)
        return read_slice(files[0]), read_slice(files[-1])

    def test_blast_states_and_short_evolution(self):
        states = {}
        for case, expected_b in ((8, 0.0), (9, 0.1)):
            with self.subTest(case=case):
                initial, final = self.run_case(case)
                states[case] = initial
                self.assertTrue(initial['density'])
                for rho in initial['density']:
                    self.assertAlmostEqual(rho / (700 * 1.67262192e-24), 1.0, delta=1.e-5)
                # Mapping and slice interpolation spread the discontinuity over
                # neighboring cells; check the resolved inner/ambient plateaus.
                self.assertAlmostEqual(max(initial['P']), 1.e-3, delta=2.e-6)
                self.assertTrue(any(abs(p / 1.e-7 - 1) < 1.e-3 for p in initial['P']))
                for name in ('Vr', 'Vt', 'Vp'):
                    self.assertTrue(all(v == 0 for v in initial[name]))
                for components in zip(initial['Br'], initial['Bt'], initial['Bp']):
                    magnitude = math.sqrt(sum(v*v for v in components))
                    self.assertAlmostEqual(magnitude, expected_b, delta=max(1.e-15, expected_b*2.e-3))
                for values in final.values():
                    self.assertTrue(all(math.isfinite(v) for v in values))
                self.assertGreater(min(final['density']), 0)
                self.assertGreater(min(final['P']), 0)
                self.assertTrue(any(abs(v) > 1 for v in final['Vr']))
                if case == 8:
                    self.assertTrue(all(v == 0 for name in ('Br', 'Bt', 'Bp') for v in final[name]))
        # Adding a uniform B must preserve the initial gas pressure. This also
        # checks background-energy subtraction and restoration for case 9.
        if len(states) == 2:
            for hydro, mhd in zip(states[8]['P'], states[9]['P']):
                self.assertAlmostEqual(hydro, mhd, delta=1.e-10)


if __name__ == '__main__':
    unittest.main()
