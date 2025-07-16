#!/bin/bash

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

# Print system version info on screen.
echo '---- System ----'
uname -a
df -h .
free -m
ulimit -a

echo '---- Compiler & env variables used ----'
echo "COMPILER=$COMPILER"
echo "MPIF90=$MPIF90"
echo "MODE=$MODE"
echo "DEFINE=$DEFINE"
echo "MPIPROCS=$MPIPROCS"
$COMPILER --version
PYTHONCMD=$(which python3)
if [ "$PYTHONCMD"x = x ]; then
  PYTHONCMD=$(which python)
fi

echo '---- Python ----'
$PYTHONCMD -V
echo -n 'numpy version: '
$PYTHONCMD -c 'import numpy; print(numpy.__version__)'
echo -n 'matplotlib version: '
$PYTHONCMD -c 'import matplotlib; print(matplotlib.__version__)'
echo -n 'sdfr version: '
$PYTHONCMD -c 'import sdfr; print(sdfr.__version__)'

