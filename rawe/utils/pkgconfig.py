# Copyright 2012-2013 Greg Horn
#
# This file is part of rawesome.
#
# rawesome is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rawesome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

import subprocess
import sys
import os

def call(args):
    # make sure pkg-config is available
    try:
        subprocess.check_call(["pkg-config", '--version'], \
                                  stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    except (subprocess.CalledProcessError, OSError):
        sys.stderr.write('`pkg-config` has not been found, please install it\n')
        sys.exit(os.EX_CONFIG)

    # make sure pkg-config can find acado and ocg2
    try:
        p = subprocess.Popen(['pkg-config']+args, \
                                 stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (output, error) = p.communicate()
        if error != '':
            raise OSError
    except (subprocess.CalledProcessError, OSError):
        raise Exception('\n'+error)

    return output.strip()

def getRpath(name):
    rpath = None
    for blah in call(['--libs',name]).split(' '):
        blah = blah.strip()
        if len(blah) > 1 and blah[:2] == '-L':
            rpath = blah[2:]
            break
    assert rpath is not None, "couldn't detect the library path of \""+name+"\" :("
    return rpath
