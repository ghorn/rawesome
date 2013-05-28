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

import mhempc_pb2
import time

def loadProtos(filename):
    sizesfilename = filename+".sizes"
    print 'reading sizes file "'+sizesfilename+'"...'
    f = open(sizesfilename, 'r')
    sizes = []
    for line in f:
        sizes.append(int(line.strip()))
    f.close()
    
    print 'reading from "'+filename+'"...'
    f = open(filename, 'r')
    raw = []
    t0 = time.time()
    for size in sizes:
        raw.append( f.read(size) )
    dt = time.time() - t0
    f.close()
    print "read %d bytes in %.3f seconds, %.3g bytes/s\n" % (sum(sizes), dt, sum(sizes)/dt)
    
    print 'deserializing protos...'
    protos = []
    t0 = time.time()
    for r in raw:
        proto = mhempc_pb2.MheMpcHorizons()
        proto.ParseFromString( r )
        protos.append( proto )
    dt = time.time() - t0
    print "desearialized %d protos in %.2f seconds, %.3g bytes/s" % (len(protos), dt, sum(sizes)/dt)

    return protos
