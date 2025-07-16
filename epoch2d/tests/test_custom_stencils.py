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

x_min = -12 * micron
x_max = -x_min
y_min = x_min
y_max = x_max

dt_multiplier = 0.99

lambda_l = 0.5 * micron
x0 = -12.0 * micron  # m
t0 = 8 * femto  # s

# derived quantities from the above
k_l = 2*np.pi/lambda_l
dx = (x_max-x_min)/nx
dy = (y_max-y_min)/ny

dt_lehe   = dt_multiplier * 1.0/np.sqrt(max(1.0/dx**2, 1.0/dy**2)) / c
dt_yee    = dt_multiplier * dx * dy / np.sqrt(dx**2 + dy**2) / c
dt_pukhov = dt_multiplier * min(dx, dy)/c

dt = dict()
dt['optimized'] = 0.9082126568805592 * dx / c
dt['optimized_symm'] = 0.8988685682513151 * dx / c
dt['optimized_xaxis'] = 0.956632159129662 * dx / c

vg = dict()
# calculate_omega.py --Y 3 --dt-multiplier 0.99 --params  0.9082126568805592
# 0.04075757835916255 -0.04142032920970152 -0.20827814817872584
# --physical-dx 1e-7 --laser-wavelength 5e-7
vg['optimized'] = 1.0490493627815458 * c

# calculate_omega.py --Y 3 --dt-multiplier 0.99 --symmetric --params
# 0.8988685682513151 0.01862292597327679 -0.04155873453935287 --physical-dx
# 1e-7 --laser-wavelength 5e-7
vg['optimized_symm'] = 1.044753207834214 * c

# calculate_omega.py --Y 3 --dt-multiplier 0.99 --params 0.956632159129662
# 0.025096871992206993 -0.017744324957063393 -0.0009692545471922645
# --physical-dx 1e-7 --laser-wavelength 5e-7
vg['optimized_xaxis'] = 1.0197513694119302 * c

vg_lehe   = c*(1.0 + 2.0*(1.0-c*dt_lehe/dx)*(k_l*dx/2.0)**2)
vg_yee    = c*np.cos(k_l*dx/2.0)/np.sqrt(1-(c*dt_yee/dx*np.sin(k_l*dx/2.0))**2)
vg_pukhov = c*np.cos(k_l*dx/2.0) \
          / np.sqrt(1-(c*dt_pukhov/dx*np.sin(k_l*dx/2.0))**2)


def xt2(sdffile, key='Electric Field/Ey'):
    t = sdffile['Header']['time']
    xaxis = sdffile[key].grid_mid.data[0]
    data = sdffile[key].data
    b = np.sum(data**2)
    if b > 0 and t > 0:
        x = np.sum(xaxis[:, np.newaxis]*data**2)/b
    else:
        x = None

    return t, x


class test_custom_stencils(SimTest):
    # solvers = ['optimized', 'optimized_symm', 'optimized_xaxis']
    solvers = ['optimized']

    @classmethod
    def setUpClass(cls):
        super(test_custom_stencils, cls).setUpClass()

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
                array = dump[key].data
                xaxis = dump[key].grid_mid.data[0]*1e6
                yaxis = dump[key].grid_mid.data[1]*1e6
                im = ax.imshow(array.T, extent=[min(xaxis), max(xaxis),
                                                min(yaxis), max(yaxis)])
                im.set_clim([-1e11, 1e11])
                t = dump['Header']['time']
                dx = (dump[key].grid_mid.data[0][1]
                      - dump[key].grid_mid.data[0][0]).item()

                vg_yee = c*np.cos(k_l*dx/2.0) \
                    / np.sqrt(1-(c*dt[solver]/dx*np.sin(k_l*dx/2.0))**2)

                x = (x0 + c*(t-t0))*1e6
                x_this = (x0 + vg[solver]*(t-t0))*1e6
                x_yee = (x0 + vg_yee*(t-t0))*1e6
                t, x_sim = xt2(dump)
                x_sim *= 1e6

                if x > 0:
                    ax.plot([x, x], [min(yaxis), 0.9*min(yaxis)],
                            color='k', linestyle='-', linewidth=2, alpha=0.25)
                    ax.plot([x, x], [0.9*max(yaxis), max(yaxis)],
                            color='k', linestyle='-', linewidth=2, alpha=0.25)

                    ax.plot([x_this, x_this], [min(yaxis), 0.9*min(yaxis)],
                            color='b', linestyle='-', linewidth=2, alpha=0.25)
                    ax.plot([x_this, x_this], [0.9*max(yaxis), max(yaxis)],
                            color='b', linestyle='-', linewidth=2, alpha=0.25)

                    ax.plot([x_sim, x_sim], [min(yaxis), 0.9*min(yaxis)],
                            color='w', linestyle='-', linewidth=2, alpha=0.25)
                    ax.plot([x_sim, x_sim], [0.9*max(yaxis), max(yaxis)],
                            color='w', linestyle='-', linewidth=2, alpha=0.25)

                    ax.plot([x_yee, x_yee], [min(yaxis), 0.9*min(yaxis)],
                            color='r', linestyle='-', linewidth=2, alpha=0.25)
                    ax.plot([x_yee, x_yee], [0.9*max(yaxis), max(yaxis)],
                            color='r', linestyle='-', linewidth=2, alpha=0.25)

                ax.set_title('{:2.1f} fs'.format(t*1e15))
                ax.set_xlabel(r'x [${\mu}\mathrm{m}$]')
                ax.axis([min(xaxis), max(xaxis), min(yaxis), max(yaxis)])

            if solver == 'optimized':
                title = 'Opt'
            elif solver == 'optimized_symm':
                title = 'Opt Symm'
            elif solver == 'optimized_xaxis':
                title = 'Opt xaxis'

            fig.suptitle(key+', c*dt/dx={:0.03f}'.format(c * dt[solver] / dx))

            fig.tight_layout()
            fig.savefig(key.replace('/', '_') + '_' + title + '.png', dpi=320)

    def test_group_velocity(self):
        tx = {}
        for solver in self.solvers:
            tx[solver] = np.array([xt2(dump)
                                   for dump in self.dumps[solver][1:]])
        print(tx)

        for solver, data in tx.items():
            vg_sim = np.polyfit(data[:, 0], data[:, 1], 1)[0]

            # For reference, right here, right now the following line prints

            # optimized       314241436.846 314497087.032 0.000812885703
            # optimized_symm  312578029.167 313209132.180 0.002014957255
            # optimized_xaxis 305472651.829 305713769.585 0.000788704270

            print('{:15} {:.3f} {:.3f} {:.12f}'.format(solver, vg_sim,
                  vg[solver], abs(vg_sim-vg[solver])/vg[solver]))

            assert np.isclose(vg_sim, vg[solver], rtol=0.003)


if __name__ == '__main__':
    unittest.main()
