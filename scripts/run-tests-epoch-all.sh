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

# THIS FILE MUST RUN WITHOUT ERROR ON EVERY COMMIT!

set -e

# chdir to epoch top level path

# The following complicated section simulates the use of "readlink -f" to find
# the source directory on platforms which don't support the "-f" flag (eg. OS X)
# It is equivalent to: cd $(dirname "$(readlink -f "$0")")/..
target="$0"
cd $(dirname "$target")
target=$(basename "$target")

# Iterate down a (possible) chain of symlinks
n=0
while [ -L "$target" ]; do
  n=$((n+1))
  if [ $n -gt 1000 ]; then
    echo "ERROR finding source directory."
    exit 1
  fi
  target=$(readlink "$target")
  cd $(dirname "$target")
  target=$(basename "$target")
done
cd ..

PYTHONCMD=$(which python3)
if [ "$PYTHONCMD"x = x ]; then
  PYTHONCMD=$(which python)
fi

# show system info
scripts/system_info.sh

# run the actual tests
$PYTHONCMD scripts/run-tests.py 1d $@
$PYTHONCMD scripts/run-tests.py 2d $@
$PYTHONCMD scripts/run-tests.py 3d $@
