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

# this makefile should be included into each makefile in the examples folder.
# its a template adding some useful targets. only in exceptions those targets
# can be changed inside the actual makefile for the example. this ensures,
# that all examples work the same, so they can be used for testing.

EPOCH?=../../bin/epoch2d
# MPRIPROCS can be set as an environmental variable, using
# $export MPRIPROCS=8
MPIPROCS?=2

default: epoch2d.dat

epoch2d.dat: input.deck $(EPOCH)
	echo . | mpirun -n $(MPIPROCS) $(EPOCH)

clean:
	rm -rf *.png *.sdf deck.status epoch1d.dat epoch2d.dat epoch3d.dat *.visit
	rm -rf __pycache__
