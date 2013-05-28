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

import kiteproto
import kite_pb2

def showAllPoints(traj,myiter,ocp,conf):
    return normalCallback(traj,myiter,ocp,conf,showAllPoints=True)

def normalCallback(traj,myiter,ocp,conf,showAllPoints=False):
    kiteProtos = []
    for k in range(0,ocp.nk):
        for nicpIdx in range(0,ocp.nicp):
            if showAllPoints:
                degIdxRange = range(ocp.deg+1)
            else:
                degIdxRange = [0]
            for degIdx in degIdxRange:
                lookup = lambda name: traj.lookup(name,timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)
                kiteProtos.append( kiteproto.toKiteProto(lookup,
                                                         lineAlpha=0.2) )
    mc = kite_pb2.MultiCarousel()
    mc.horizon.extend(list(kiteProtos))

    mc.messages.append("w0: "+str(traj.lookup('w0')))
    mc.messages.append("iter: "+str(myiter))
    mc.messages.append("endTime: "+str(traj.lookup('endTime')))
    mc.messages.append("average electrical power: "+str(traj.lookup('electrical_energy',timestep=-1)/traj.lookup('endTime'))+" W")
    mc.messages.append("average mechanical power: "+str(traj.lookup('mechanical_energy',timestep=-1)/traj.lookup('endTime'))+" W")
    
    return mc.SerializeToString()
