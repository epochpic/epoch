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

import numpy as np
import sdfr
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import unittest
import platform
from . import SimTest


def showdatafields(sdffile):
    data = sdfr.read(sdffile, dict=True)
    print('Data fields present in "' + sdffile + '":')
    print(data.keys())


def plotdump(sdffile, key, ax):
    data = sdfr.read(sdffile, dict=True)
    array = data[key].data
    ax.plot(array)
    ax.set_title(r'{:2.1f} ms'.format(data['Header']['time']*1e3))


def plotdump2d(sdffile, key, ax):
    data = sdfr.read(sdffile, dict=True)
    if len(data[key].dims) == 3:
        array = data[key].data[:, :, 0]
    else:
        array = data[key].data[:, :]
    ax.imshow(array)
    ax.set_title(r'{:2.1f} ms'.format(data['Header']['time']*1e3))


def plotevolution(key):
    print('plotting: ' + key)
    sdffiles = ['{:04d}.sdf'.format(i) for i in range(15)]
    fig, axarr = plt.subplots(3, 5, figsize=(20, 11))
    for (sdffile, ax) in zip(sdffiles, np.ravel(axarr)):
        plotdump(sdffile, key, ax)
    fig.suptitle(key)
    fig.savefig(key.replace('/', '_') + '.png', dpi=160)


def plot2devolution(key):
    print('plotting: ' + key)
    sdffiles = ['{:04d}.sdf'.format(i) for i in range(15)]
    fig, axarr = plt.subplots(3, 5, figsize=(20, 11))
    for (sdffile, ax) in zip(sdffiles, np.ravel(axarr)):
        plotdump2d(sdffile, key, ax)
    fig.suptitle(key)
    fig.savefig(key.replace('/', '_') + '.png', dpi=160)


def createplots():
    if platform.system() == 'Darwin':
        print('macosx backend')
        plt.switch_backend('macosx')
    # showdatafields('0000.sdf')
    plotevolution('Electric Field/Ex')
    plotevolution('Derived/Number_Density/electrons')
    plotevolution('Derived/Number_Density/protons')
    plotevolution('Derived/Charge_Density')
    plot2devolution('dist_fn/x_px/electrons')
    plot2devolution('dist_fn/x_px/protons')


class test_landau(SimTest):

    def test_creatplots(self):
        createplots()


if __name__ == '__main__':
    unittest.main()
