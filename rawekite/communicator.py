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

import zmq
import kiteproto

class Communicator(object):
    def __init__(self):
        self.context   = zmq.Context(1)
        self.publisher = self.context.socket(zmq.PUB)
        self.publisher.bind("tcp://*:5563")
#        self.fOutputs = fOutputs
#        self.outputNames = outputNames

    def sendKite(self,x,u,p,outs,conf,otherMessages=[]):
        pb = kiteproto.toKiteProto(dict(x.items()+u.items()+p.items()+outs.items()).__getitem__,
                                   lineAlpha=0.2)

        for name in outs:
            v = outs[name]
            if type(v) is float:
                pb.messages.append(name+": "+str(v))
        if len(otherMessages)>0:
            pb.messages.append("-------------------------")
            for om in otherMessages:
                pb.messages.append(om)
        self.publisher.send_multipart(["carousel", pb.SerializeToString()])
        
    def close(self):
        self.publisher.close()
        self.context.term()
