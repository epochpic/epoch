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
Shared helpers for the epoch3d custom-laser tests.

Files are raw float64, tr1 fastest, then tr2, time slowest, which is
what tofile() writes for a C-ordered (n_t, n_tr2, n_tr1) array.
'''

import glob
import os
import subprocess

import numpy as np


def write_spatiotemporal(path, array_t_tr2_tr1):
    '''Write an (n_t, n_tr2, n_tr1) array.'''
    np.asarray(array_t_tr2_tr1, dtype=np.float64).tofile(str(path))


def write_static_plane(path, array_tr2_tr1):
    '''Write an (n_tr2, n_tr1) array.'''
    np.asarray(array_tr2_tr1, dtype=np.float64).tofile(str(path))


def gaussian_2d(tr1, tr2, tr1_0, tr2_0, w1, w2, amp=1.0):
    '''2D Gaussian on the tr1 x tr2 grid, shape (len(tr2), len(tr1)).'''
    t1, t2 = np.meshgrid(tr1, tr2)
    return amp * np.exp(-((t1 - tr1_0) ** 2) / w1 ** 2
                         - ((t2 - tr2_0) ** 2) / w2 ** 2)


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
    patterns = ('*.sdf', '*.visit', '*.png', 'deck.status', 'epoch3d.dat')
    for pattern in patterns:
        for f in glob.glob(os.path.join(str(run_dir), pattern)):
            os.remove(f)
