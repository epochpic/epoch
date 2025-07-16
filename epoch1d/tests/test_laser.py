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
from tests import SimTest


def showdatafields(sdffile):
    data = sdfr.read(sdffile, dict=True)
    print('Data fields present in "' + sdffile + '":')
    print(data.keys())


def plotdump(sdffile, key, ax):
    data = sdfr.read(sdffile, dict=True)
    array = data[key].data
    ax.plot(array)
    ax.set_title('{:2.1f} fs'.format(data['Header']['time']*1e15))


def plotevolution(key):
    sdffiles = ['{:04d}.sdf'.format(i) for i in range(8)]
    fig, axarr = plt.subplots(2, 4, figsize=(16, 9))
    for (sdffile, ax) in zip(sdffiles, np.ravel(axarr)):
        plotdump(sdffile, key, ax)
    fig.suptitle(key)
    fig.savefig(key.replace('/', '_') + '.png', dpi=160)


def createplots():
    if platform.system() == 'Darwin':
        print('macosx backend')
        plt.switch_backend('macosx')
    # showdatafields('0000.sdf')
    plotevolution('Electric Field/Ey')


class test_laser(SimTest):

    def totaleyassert(self, dump, val):
        data = sdfr.read(dump, dict=True)
        ey = data['Electric Field/Ey'].data**2
        eysum = np.sum(ey)
        test = np.isclose(eysum, val)
        msg = '{:} != {:} within the given tolerances.'.format(eysum, val)
        self.assertTrue(test, msg=msg)

    def test_createplots(self):
        createplots()

    def test_Eey0000(self):
        self.totaleyassert('0000.sdf', 0.0)

    def test_Eey0001(self):
        self.totaleyassert('0001.sdf', 1.38636e+23)

    def test_Eey0003(self):
        self.totaleyassert('0003.sdf', 1.40618e+23)

    def test_Eey0007(self):
        self.totaleyassert('0007.sdf', 6.90067e+17)


if __name__ == '__main__':
    unittest.main()
