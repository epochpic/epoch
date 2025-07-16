#!/usr/bin/env python

# Copyright (C) 2009-2019 University of Warwick
# Stephan Kuschel <Stephan.Kuschel@gmail.com>
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

import os.path as osp
import unittest
import platform

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sdfr

from . import SimTest


micron = 1e-6
femto = 1e-15
c = 2.99792458e8  # m/s

# check that these correspond to input deck!

nx = 240
ny = nx//3
nz = nx//3

x_min = -12 * micron
x_max = -x_min
y_min = x_min
y_max = x_max
z_min = x_min
z_max = x_max

dt_multiplier = 0.95

lambda_l = 0.5 * micron
x0 = -12.0 * micron  # m
t0 = 8 * femto  # s

# derived quantities from the above
k_l = 2*np.pi/lambda_l
dx = (x_max-x_min)/nx
dy = (y_max-y_min)/ny
dz = (z_max-z_min)/nz

dt_lehe   = dt_multiplier \
          * 1.0/np.sqrt(max(1.0/dx**2, 1.0/dy**2 + 1.0/dz**2)) / c
dt_yee    = dt_multiplier * dx * dy * dz \
          / np.sqrt((dx*dy)**2 + (dy*dz)**2 + (dz*dx)**2) / c
dt_pukhov = dt_multiplier * min(dx, dy, dz)/c
dt_cowan  = dt_pukhov

vg_lehe   = c*(1.0 + 2.0*(1.0-c*dt_lehe/dx)*(k_l*dx/2.0)**2)
vg_yee    = c*np.cos(k_l*dx/2.0)/np.sqrt(1-(c*dt_yee/dx*np.sin(k_l*dx/2.0))**2)
vg_pukhov = c*np.cos(k_l*dx/2.0) \
          / np.sqrt(1-(c*dt_pukhov/dx*np.sin(k_l*dx/2.0))**2)
vg_cowan  = c*np.cos(k_l*dx/2.0) \
          / np.sqrt(1-(c*dt_cowan/dx*np.sin(k_l*dx/2.0))**2)


def xt3(sdffile, key='Electric Field/Ey'):
    t = sdffile['Header']['time']
    xaxis = sdffile[key].grid_mid.data[0]
    data = sdffile[key].data
    b = np.sum(data**2)
    if b > 0 and t > 0:
        x = np.sum(xaxis[:, np.newaxis, np.newaxis]*data**2)/b
    else:
        x = None

    return t, x


class test_maxwell_solvers(SimTest):
    # solvers = ['yee', 'lehe_x', 'pukhov', 'cowan']
    solvers = ['lehe_x', 'cowan']

    @classmethod
    def setUpClass(cls):
        super(test_maxwell_solvers, cls).setUpClass()

        dumps = cls.dumps = {}
        for solver in cls.solvers:
            l = dumps.setdefault(solver, [])
            for dump in [osp.join(solver, '{:04d}.sdf'.format(i))
                         for i in range(4)]:
                l.append(sdfr.read(dump, dict=True))

    def test_createplot(self):
        if platform.system() == 'Darwin':
            print('macosx backend')
            plt.switch_backend('macosx')

        for solver in self.solvers:
            key = 'Electric Field/Ey'
            fig, axarr = plt.subplots(1, 3, figsize=(10, 4))
            axarr[0].set_ylabel(r'y [${\mu}\mathrm{m}$]')
            dumps = self.dumps[solver][1:]
            for (dump, ax) in zip(dumps, np.ravel(axarr)):
                array = dump[key].data[:, nz//2, :]
                xaxis = dump[key].grid_mid.data[0]*1e6
                yaxis = dump[key].grid_mid.data[1]*1e6
                im = ax.imshow(array.T, extent=[min(xaxis), max(xaxis),
                                                min(yaxis), max(yaxis)])
                im.set_clim([-1e11, 1e11])
                t = dump['Header']['time']

                x = (x0 + c*(t-t0))*1e6
                x_lehe = (x0 + vg_lehe*(t-t0))*1e6
                x_yee = (x0 + vg_yee*(t-t0))*1e6

                if x > 0:
                    ax.plot([x, x], [min(yaxis), 0.9*min(yaxis)], color='k',
                            linestyle='-', linewidth=2, alpha=0.25)
                    ax.plot([x, x], [0.9*max(yaxis), max(yaxis)], color='k',
                            linestyle='-', linewidth=2, alpha=0.25)

                    ax.plot([x_lehe, x_lehe], [min(yaxis), 0.9*min(yaxis)],
                            color='b', linestyle='-', linewidth=2, alpha=0.25)
                    ax.plot([x_lehe, x_lehe], [0.9*max(yaxis), max(yaxis)],
                            color='b', linestyle='-', linewidth=2, alpha=0.25)

                    ax.plot([x_yee, x_yee], [min(yaxis), 0.9*min(yaxis)],
                            color='r', linestyle='-', linewidth=2, alpha=0.25)
                    ax.plot([x_yee, x_yee], [0.9*max(yaxis), max(yaxis)],
                            color='r', linestyle='-', linewidth=2, alpha=0.25)

                ax.set_title('{:2.1f} fs'.format(dump['Header']['time']*1e15))

            if solver == 'lehe_x':
                dt = c * dt_lehe / dx
                title = 'Lehe'
            elif solver == 'pukhov':
                dt = c * dt_pukhov / dx
                title = 'Pukhov'
            elif solver == 'yee':
                dt = c * dt_yee / dx
                title = 'Yee'
            elif solver == 'cowan':
                dt = c * dt_cowan / dx
                title = 'Cowan'

            fig.suptitle(key+', c*dt/dx={:0.03f}'.format(dt))

            fig.tight_layout()
            fig.savefig(key.replace('/', '_') + '_' + title + '.png', dpi=320)

    def test_group_velocity(self):
        tx = {}
        for solver in self.solvers:
            tx[solver] = np.array([xt3(dump)
                                   for dump in self.dumps[solver][1:]])
        print(tx)

        vg = dict(lehe_x=vg_lehe, pukhov=vg_pukhov, yee=vg_yee, cowan=vg_cowan)

        for solver, data in tx.items():
            vg_sim = np.polyfit(data[:, 0], data[:, 1], 1)[0]

            # For reference, right here, right now the following line prints

            # pukhov 291961761.344 292363351.796 0.001373600521
            # yee    279231545.324 281017561.535 0.006355532378
            # cowan  291962719.038 292363351.796 0.001370324821
            # lehe_x 311952693.446 311627789.852 0.001042601480

            print('{:6} {:.3f} {:.3f} {:.12f}'.format(solver, vg_sim,
                  vg[solver], abs(vg_sim-vg[solver])/vg[solver]))

            assert np.isclose(vg_sim, vg[solver], rtol=0.007)


if __name__ == '__main__':
    unittest.main()
