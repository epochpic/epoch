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
Shared helpers for the epoch2d custom-laser tests.

Files are raw float64, transverse axis fastest and time slowest, which
is what tofile() writes for a C-ordered (n_t, n_transverse) array.
'''

import glob
import os
import subprocess

import numpy as np


def write_spatiotemporal(path, array_t_tr):
    '''Write an (n_t, n_transverse) array.'''
    np.asarray(array_t_tr, dtype=np.float64).tofile(str(path))


def write_static_line(path, array_tr):
    '''Write an (n_transverse,) array.'''
    np.asarray(array_tr, dtype=np.float64).tofile(str(path))


def gaussian_1d(tr, tr_0, w, amp=1.0):
    '''Gaussian centred on tr_0 with 1/e half-width w.'''
    return amp * np.exp(-((tr - tr_0) ** 2) / w ** 2)


def pick_active_snapshot(dumps, key, t_start, t_end):
    '''The dump in [t_start, t_end] with the largest peak |key|. A fixed
    time could land on a carrier zero-crossing, where field ratios are
    meaningless.'''
    candidates = [d for d in dumps if t_start <= d['Header']['time'] <= t_end]
    scores = [np.abs(d[key].data).max() for d in candidates]
    return candidates[int(np.argmax(scores))]


def run_epoch(run_dir, epoch_bin, nprocs=2, timeout=180):
    '''Run epoch directly, as makefile.inc does, and return (exit code,
    output). For tests that expect it to fail.'''
    proc = subprocess.run(
        ['mpirun', '-n', str(nprocs), str(epoch_bin)],
        cwd=str(run_dir), input='.\n', capture_output=True, text=True,
        timeout=timeout)
    return proc.returncode, proc.stdout + '\n' + proc.stderr


def run_make(sim_dir):
    '''Run 'make clean' then 'make' in sim_dir; return make's exit code.
    The clean is needed because the make target depends only on
    input.deck, so regenerated fixtures alone would not force a rerun.'''
    cwd = os.getcwd()
    try:
        os.chdir(str(sim_dir))
        subprocess.call('make clean', shell=True)
        return subprocess.call('make', shell=True)
    finally:
        os.chdir(cwd)


def run_make_clean(sim_dir):
    '''Run 'make clean' in sim_dir unless NOCLEAN is set.'''
    if 'NOCLEAN' in os.environ:
        return
    cwd = os.getcwd()
    try:
        os.chdir(str(sim_dir))
        subprocess.call('make clean', shell=True)
    finally:
        os.chdir(cwd)


def remove_run_artifacts(run_dir):
    '''Delete the outputs 'make clean' would, for runs started with
    run_epoch.'''
    patterns = ('*.sdf', '*.visit', '*.png', 'deck.status', 'epoch2d.dat')
    for pattern in patterns:
        for f in glob.glob(os.path.join(str(run_dir), pattern)):
            os.remove(f)
