#!/usr/bin/env python

# Copyright (C) 2009-2019 University of Warwick
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Regression tests for custom laser profile injection in epoch2d (see
DOCUMENTATION_LASER_INJECTION.pdf). As epoch3d's tests, but with one
transverse axis and without the epoch3d-only slab and alias cases.

Each class writes its own profile/phase fixtures, runs epoch2d in
tests/custom_laser/<case>/ and checks one property of the result.
Field comparisons use common.pick_active_snapshot, since a fixed time
can land on a carrier zero-crossing. test_custom_laser_size_abort
expects epoch2d to abort, so it runs it directly via common.run_epoch.
'''

import os
import unittest

import numpy as np
import sdfr

from .custom_laser import common as cl

micron = 1e-6
femto = 1e-15

HERE = os.path.dirname(__file__)
CUSTOM_LASER_DIR = os.path.join(HERE, 'custom_laser')
EPOCH_BIN = os.path.join(HERE, '..', 'bin', 'epoch2d')


def dump_series(simdir):
    '''All *.sdf dumps in simdir, in time order.'''
    files = sorted(f for f in os.listdir(simdir) if f.endswith('.sdf'))
    return [sdfr.read(os.path.join(simdir, f), dict=True) for f in files]


class _EpochRanCase(unittest.TestCase):
    '''Fails each test up front if setUpClass's epoch2d run failed, as
    SimTest does for single-level test directories.'''

    def setUp(self):
        if self.epochexitcode:
            self.fail('running epoch2d errored (exitcode {})'
                      .format(self.epochexitcode))


class test_custom_laser_y_min_boundary(_EpochRanCase):
    '''Spatiotemporal laser on y_min: the transverse axis is x (the
    deck still names its bounds y_min/y_max, a documented 2D quirk) and
    pol_angle = 0 should drive Ez. A Gaussian spot drifting in x is
    tracked against the file.'''

    N_TR, N_T = 10, 6
    SIMDIR = os.path.join(CUSTOM_LASER_DIR, 'y_min_boundary')

    @classmethod
    def setUpClass(cls):
        x_file = np.linspace(-4, 4, cls.N_TR) * micron
        x0_file = np.linspace(-3, 3, cls.N_T) * micron
        w = 0.6 * micron
        arr = np.zeros((cls.N_T, cls.N_TR))
        for it in range(cls.N_T):
            arr[it] = cl.gaussian_1d(x_file, x0_file[it], w)
        cl.write_spatiotemporal(
            os.path.join(cls.SIMDIR, 'profile_y_min.dat'), arr)

        cls.epochexitcode = cl.run_make(cls.SIMDIR)
        if not cls.epochexitcode:
            cls.dumps = dump_series(cls.SIMDIR)

    @classmethod
    def tearDownClass(cls):
        cl.run_make_clean(cls.SIMDIR)
        fixture = os.path.join(cls.SIMDIR, 'profile_y_min.dat')
        if os.path.exists(fixture):
            os.remove(fixture)

    def test_ez_channel_dominates_ex(self):
        '''pol_angle = 0 on a y_min laser should drive Ez, not Ex.'''
        ex_max = max(np.abs(d['Electric Field/Ex'].data).max()
                      for d in self.dumps)
        ez_max = max(np.abs(d['Electric Field/Ez'].data).max()
                      for d in self.dumps)
        self.assertGreater(ez_max, 10.0 * ex_max)

    def test_peak_position_tracks_encoded_drift(self):
        '''The Ez peak's x-position should follow the drift encoded in
        the file.'''
        t_file = np.linspace(0, 8, self.N_T) * femto
        x0_file = np.linspace(-3, 3, self.N_T) * micron
        x_grid = self.dumps[0]['Electric Field/Ez'].grid_mid.data[0]

        peaks, times = [], []
        for d in self.dumps:
            t = d['Header']['time']
            if t < 1.0 * femto or t > 8.0 * femto:
                continue
            col = d['Electric Field/Ez'].data[:, 0]
            ix = int(np.argmax(np.abs(col)))
            peaks.append(x_grid[ix])
            times.append(t)

        encoded = np.interp(times, t_file, x0_file)
        corr = np.corrcoef(peaks, encoded)[0, 1]
        self.assertGreater(corr, 0.9)


class test_custom_laser_two_laser_independence(_EpochRanCase):
    '''Two x_min lasers (pol = 0 and 90) with different files must not
    share loaded data: Ey and Ez should peak where their own files put
    them, on opposite sides of y = 0.'''

    N_TR, N_T = 10, 5
    T_START, T_END = 0.0, 6 * femto
    SIMDIR = os.path.join(CUSTOM_LASER_DIR, 'two_laser_independence')

    @classmethod
    def setUpClass(cls):
        y_file = np.linspace(-4, 4, cls.N_TR) * micron
        cls.peak_a_y = y_file[7]
        cls.peak_b_y = y_file[2]

        a_line = cl.gaussian_1d(y_file, cls.peak_a_y, 0.6 * micron,
                                 amp=1.0)
        cl.write_spatiotemporal(os.path.join(cls.SIMDIR, 'profile_A.dat'),
                                 np.tile(a_line, (cls.N_T, 1)))

        b_line = cl.gaussian_1d(y_file, cls.peak_b_y, 1.8 * micron,
                                 amp=0.4)
        cl.write_spatiotemporal(os.path.join(cls.SIMDIR, 'profile_B.dat'),
                                 np.tile(b_line, (cls.N_T, 1)))

        cls.epochexitcode = cl.run_make(cls.SIMDIR)
        if not cls.epochexitcode:
            dumps = dump_series(cls.SIMDIR)
            cls.snap = cl.pick_active_snapshot(
                dumps, 'Electric Field/Ey', cls.T_START, cls.T_END)

    @classmethod
    def tearDownClass(cls):
        cl.run_make_clean(cls.SIMDIR)
        for name in ('profile_A.dat', 'profile_B.dat'):
            fixture = os.path.join(cls.SIMDIR, name)
            if os.path.exists(fixture):
                os.remove(fixture)

    def test_channels_independent(self):
        y_grid = self.snap['Electric Field/Ey'].grid_mid.data[1]
        ey_line = self.snap['Electric Field/Ey'].data[0, :]
        ez_line = self.snap['Electric Field/Ez'].data[0, :]

        iy_ey = int(np.argmax(np.abs(ey_line)))
        iy_ez = int(np.argmax(np.abs(ez_line)))

        dy = y_grid[1] - y_grid[0]
        self.assertLess(abs(y_grid[iy_ey] - self.peak_a_y), 1.5 * dy)
        self.assertLess(abs(y_grid[iy_ez] - self.peak_b_y), 1.5 * dy)
        # Shared (aliased) storage would put both peaks on the same side.
        self.assertLess(y_grid[iy_ey] * y_grid[iy_ez], 0.0)


class test_custom_laser_interpolation_accuracy(_EpochRanCase):
    '''Checks the injected field against a Python copy of
    custom_laser_profile. The file is separable, g(t) * h(y), so field
    ratios between points equal ratios of h whatever the carrier
    phase.'''

    N_TR, N_T = 8, 4
    TR_MIN, TR_MAX = -4 * micron, 4 * micron
    T_START, T_END = 0.0, 6 * femto
    Y0, WY = 1.0 * micron, 2.0 * micron
    SIMDIR = os.path.join(CUSTOM_LASER_DIR, 'interpolation_accuracy')

    @classmethod
    def setUpClass(cls):
        y_file = np.linspace(cls.TR_MIN, cls.TR_MAX, cls.N_TR)
        t_file = np.linspace(cls.T_START, cls.T_END, cls.N_T)

        def h_shape(y):
            return np.exp(-((y - cls.Y0) ** 2) / cls.WY ** 2)

        def g_time(t_norm):
            return 0.4 + 0.6 * t_norm

        arr = np.zeros((cls.N_T, cls.N_TR))
        for it, t in enumerate(t_file):
            arr[it] = g_time(it / (cls.N_T - 1)) * h_shape(y_file)
        cl.write_spatiotemporal(os.path.join(cls.SIMDIR, 'profile.dat'),
                                 arr)

        cls.y_file, cls.arr = y_file, arr

        cls.epochexitcode = cl.run_make(cls.SIMDIR)
        if not cls.epochexitcode:
            dumps = dump_series(cls.SIMDIR)
            cls.snap = cl.pick_active_snapshot(
                dumps, 'Electric Field/Ey', cls.T_START, cls.T_END)

    @classmethod
    def tearDownClass(cls):
        cl.run_make_clean(cls.SIMDIR)
        fixture = os.path.join(cls.SIMDIR, 'profile.dat')
        if os.path.exists(fixture):
            os.remove(fixture)

    def _sample_py(self, pos, time):
        '''Python copy of custom_laser_profile (0-indexed).'''
        if pos < self.TR_MIN or pos > self.TR_MAX:
            return 0.0
        if time < self.T_START or time > self.T_END:
            return 0.0
        d1 = (self.TR_MAX - self.TR_MIN) / (self.N_TR - 1)
        dt = (self.T_END - self.T_START) / (self.N_T - 1)
        i1 = int((pos - self.TR_MIN) / d1)
        it = int((time - self.T_START) / dt)
        i1 = max(0, min(i1, self.N_TR - 2))
        it = max(0, min(it, self.N_T - 2))
        p10 = self.TR_MIN + i1 * d1
        t0 = self.T_START + it * dt
        u = (pos - p10) / d1
        v = (time - t0) / dt
        m = self.arr
        f0 = (1 - u) * m[it, i1] + u * m[it, i1 + 1]
        f1 = (1 - u) * m[it + 1, i1] + u * m[it + 1, i1 + 1]
        return (1 - v) * f0 + v * f1

    def test_interpolation_matches_reference(self):
        t_snap = self.snap['Header']['time']
        y_grid = self.snap['Electric Field/Ey'].grid_mid.data[1]
        line = self.snap['Electric Field/Ey'].data[0, :]

        points = {
            'node3': float(self.y_file[3]),
            'node5': float(self.y_file[5]),
            'offnode1': float(self.y_file[2])
                + 0.3 * (self.y_file[3] - self.y_file[2]),
            'offnode2': float(self.y_file[6])
                - 0.4 * (self.y_file[6] - self.y_file[5]),
            # Not exactly TR_MAX: its nearest grid point can fall
            # just outside the file grid, where the field is zero.
            'near_edge': self.TR_MAX - 0.3 * micron,
        }
        rows = {}
        for name, y_pt in points.items():
            iy = int(np.argmin(np.abs(y_grid - y_pt)))
            y_actual = float(y_grid[iy])
            measured = float(line[iy])
            predicted = self._sample_py(y_actual, t_snap)
            rows[name] = (predicted, measured)

        ref_pred, ref_meas = rows['offnode1']
        for name, (pred, meas) in rows.items():
            pred_ratio = pred / ref_pred
            meas_ratio = meas / ref_meas
            rel_diff = abs(pred_ratio - meas_ratio) / max(abs(pred_ratio),
                                                            1e-9)
            self.assertLess(rel_diff, 0.2,
                             msg='{}: predicted/measured ratio mismatch '
                                 '({:.1%})'.format(name, rel_diff))

    def test_outside_declared_grid_is_near_zero(self):
        '''The field outside the declared file grid should be
        negligible (custom_laser_profile returns 0 there).'''
        y_grid = self.snap['Electric Field/Ey'].grid_mid.data[1]
        line = self.snap['Electric Field/Ey'].data[0, :]

        i_node5 = int(np.argmin(np.abs(y_grid - self.y_file[5])))
        i_outside = int(np.argmin(np.abs(y_grid
                                          - (self.TR_MAX + 1.0 * micron))))
        node5_val = abs(float(line[i_node5]))
        outside_val = abs(float(line[i_outside]))
        self.assertLess(outside_val, 0.05 * node5_val)


class test_custom_laser_edge_behaviour(_EpochRanCase):
    '''Outside the declared file grid, the static path clamps to the
    edge value and the spatiotemporal path returns zero. The file covers
    +-2 um of a +-6 um domain: static on y_min (Ez),
    spatiotemporal on x_min (Ey).'''

    N = 6
    FMIN, FMAX = -2 * micron, 2 * micron
    T_START, T_END = 0.0, 5 * femto
    SIMDIR = os.path.join(CUSTOM_LASER_DIR, 'edge_behaviour')

    @classmethod
    def setUpClass(cls):
        c1 = np.linspace(cls.FMIN, cls.FMAX, cls.N)
        w = 3.0 * micron
        line = cl.gaussian_1d(c1, 0.0, w, amp=1.0)
        cl.write_static_line(
            os.path.join(cls.SIMDIR, 'static', 'spatial_profile.dat'),
            line)
        cl.write_spatiotemporal(
            os.path.join(cls.SIMDIR, 'dynamic', 'profile.dat'),
            np.tile(line, (2, 1)))

        cls.epochexitcode = cl.run_make(cls.SIMDIR)
        if not cls.epochexitcode:
            static_dumps = dump_series(os.path.join(cls.SIMDIR, 'static'))
            dynamic_dumps = dump_series(
                os.path.join(cls.SIMDIR, 'dynamic'))
            # The static laser is always on; the dynamic one only
            # within its time window.
            cls.snap_static = cl.pick_active_snapshot(
                static_dumps, 'Electric Field/Ez',
                -float('inf'), float('inf'))
            cls.snap_dynamic = cl.pick_active_snapshot(
                dynamic_dumps, 'Electric Field/Ey',
                cls.T_START, cls.T_END)

    @classmethod
    def tearDownClass(cls):
        cl.run_make_clean(cls.SIMDIR)
        for fixture in (
                os.path.join(cls.SIMDIR, 'static', 'spatial_profile.dat'),
                os.path.join(cls.SIMDIR, 'dynamic', 'profile.dat')):
            if os.path.exists(fixture):
                os.remove(fixture)

    def _value_at(self, grid, line, coord):
        i = int(np.argmin(np.abs(grid - coord)))
        return abs(float(line[i]))

    def test_static_path_clamps_to_edge_value(self):
        x_grid = self.snap_static['Electric Field/Ez'].grid_mid.data[0]
        line = self.snap_static['Electric Field/Ez'].data[:, 0]

        edge = self._value_at(x_grid, line, 2 * micron)
        beyond = self._value_at(x_grid, line, 4 * micron)
        self.assertGreater(beyond / edge, 0.5)

    def test_spatiotemporal_path_zero_pads(self):
        y_grid = self.snap_dynamic['Electric Field/Ey'].grid_mid.data[1]
        line = self.snap_dynamic['Electric Field/Ey'].data[0, :]

        edge = self._value_at(y_grid, line, 2 * micron)
        beyond = self._value_at(y_grid, line, 4 * micron)
        self.assertLess(beyond / edge, 0.05)


class test_custom_laser_phase_override(_EpochRanCase):
    '''With use_phase_from_file = T, a deck 'phase' expression must not
    override the file phase on the static path. A (tilt + file) must
    equal B (file only); C (tilt, no file) is a positive control that
    must differ.'''

    N = 6
    FMIN, FMAX = -4 * micron, 4 * micron
    CASES = ('run_A_tilt_and_file', 'run_B_file_only', 'run_C_tilt_no_file')
    SIMDIR = os.path.join(CUSTOM_LASER_DIR, 'phase_override')

    @classmethod
    def setUpClass(cls):
        flat_amp = np.ones(cls.N)
        flat_phase = np.zeros(cls.N)
        for name in cls.CASES:
            cl.write_static_line(
                os.path.join(cls.SIMDIR, name, 'spatial_profile.dat'),
                flat_amp)
            cl.write_static_line(
                os.path.join(cls.SIMDIR, name, 'phase_profile.dat'),
                flat_phase)

        cls.epochexitcode = cl.run_make(cls.SIMDIR)
        if cls.epochexitcode:
            return
        cls.profiles = {}
        for name in cls.CASES:
            last = dump_series(os.path.join(cls.SIMDIR, name))[-1]
            # Copy: sdfr arrays are freed along with their dump.
            cls.profiles[name] = last['Electric Field/Ey'].data[0, :].copy()

    @classmethod
    def tearDownClass(cls):
        cl.run_make_clean(cls.SIMDIR)
        for name in cls.CASES:
            for fname in ('spatial_profile.dat', 'phase_profile.dat'):
                fixture = os.path.join(cls.SIMDIR, name, fname)
                if os.path.exists(fixture):
                    os.remove(fixture)

    def test_file_phase_cannot_be_overridden_by_deck_expression(self):
        a = self.profiles['run_A_tilt_and_file']
        b = self.profiles['run_B_file_only']
        c = self.profiles['run_C_tilt_no_file']

        max_a = np.max(np.abs(a))
        rel_ab = np.max(np.abs(a - b)) / (max_a + 1e-300)
        rel_ac = np.max(np.abs(a - c)) / (max_a + 1e-300)

        self.assertLess(rel_ab, 1e-6,
                         msg='deck phase expression changed the result '
                             'even though use_phase_from_file = T')
        # The tilt must have a visible effect, or A == B proves nothing.
        self.assertGreater(rel_ac, 0.5)


class test_custom_laser_size_abort(unittest.TestCase):
    '''A file whose size does not match n_transverse_points *
    n_t_points * 8 bytes must abort, reporting both byte counts. Each
    count is checked, plus a correctly-sized control.'''

    N_TR, N_T = 8, 5
    SIMDIR = os.path.join(CUSTOM_LASER_DIR, 'size_abort')
    CASES = {
        'run_correct':   (N_TR,     N_T),
        'run_wrong_ny':  (N_TR + 1, N_T),
        'run_wrong_nt':  (N_TR,     N_T + 1),
    }

    @classmethod
    def setUpClass(cls):
        arr = np.random.default_rng(0).normal(size=(cls.N_T, cls.N_TR))
        cls.true_bytes = cls.N_TR * cls.N_T * 8

        cls.logs = {}
        cls.exitcodes = {}
        for name in cls.CASES:
            rundir = os.path.join(cls.SIMDIR, name)
            cl.write_spatiotemporal(os.path.join(rundir, 'profile.dat'),
                                     arr)
            rc, log = cl.run_epoch(rundir, EPOCH_BIN, nprocs=2, timeout=60)
            cls.exitcodes[name] = rc
            cls.logs[name] = log

    @classmethod
    def tearDownClass(cls):
        for name in cls.CASES:
            rundir = os.path.join(cls.SIMDIR, name)
            fixture = os.path.join(rundir, 'profile.dat')
            if os.path.exists(fixture):
                os.remove(fixture)
            if 'NOCLEAN' not in os.environ:
                cl.remove_run_artifacts(rundir)

    def test_correctly_sized_file_succeeds(self):
        self.assertEqual(self.exitcodes['run_correct'], 0)

    def test_mismatched_file_aborts_with_byte_counts(self):
        for name, (ntr, nt) in self.CASES.items():
            if name == 'run_correct':
                continue
            with self.subTest(case=name):
                declared_bytes = ntr * nt * 8
                self.assertNotEqual(self.exitcodes[name], 0)
                log = self.logs[name]
                self.assertIn(str(self.true_bytes), log)
                self.assertIn(str(declared_bytes), log)


if __name__ == '__main__':
    unittest.main()
