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
import select

def call(args, cwd = '.'):
    p = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = cwd)

    msgs = []

    while True:
        reads = [p.stdout.fileno(), p.stderr.fileno()]
        ret = select.select(reads, [], [])

        for fd in ret[0]:
            if fd == p.stdout.fileno():
                read = p.stdout.readline()
#                 sys.stdout.write(read)
                msgs.append(read)
            if fd == p.stderr.fileno():
                read = p.stderr.readline()
#                 sys.stderr.write(read)
                msgs.append(read)

        ret = p.poll()
        if ret != None:
            return (ret, (''.join(msgs)).strip())
