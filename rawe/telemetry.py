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
from multiprocessing import Manager, Process
import numpy

from rawe.collocation import trajectory

def startTelemetry(ocp, callbacks=[],
                   printBoundViolation=False,printConstraintViolation=False,
                   url="tcp://*:5563"):
    xOptQueue = Manager().Queue()
    Sender(ocp, xOptQueue, callbacks, url).start()

    class MyCallback:
        def __call__(self,f,*args):
            xOpt = numpy.array(f.input('x'))

            if printBoundViolation:
                lbx = numpy.array(ocp.solver.input('lbx'))
                ubx = numpy.array(ocp.solver.input('ubx'))
                ocp._bounds.printBoundsFeedback(xOpt,lbx,ubx,reportThreshold=0)

            if printConstraintViolation:
                lbg = numpy.array(ocp.solver.input('lbg'))
                ubg = numpy.array(ocp.solver.input('ubg'))
                ocp._gfcn.setInput(xOpt,0)
                ocp._gfcn.evaluate()
                g = ocp._gfcn.output()
                ocp._constraints.printViolations(g,lbg,ubg,reportThreshold=0)

            xOptQueue.put(xOpt)

    return MyCallback()


class Sender(Process):
    def __init__(self, ocp, xOptQueue, callbackFunsWithChannels, url):
        Process.__init__(self)
        self.daemon = True

        self.ocp = ocp
        self.xOptQueue = xOptQueue
        self.url = url
        msg = "callbacks must be a list of (callback (function), channel (string)) tuples"
        if not isinstance(callbackFunsWithChannels, list):
            callbackFunsWithChannels = [callbackFunsWithChannels]
        for cfc in callbackFunsWithChannels:
            assert isinstance(cfc, tuple), msg
            assert isinstance(cfc[1], basestring), msg
        self.callbackFunsWithChannels = callbackFunsWithChannels

    def run(self):
        myiter = 0
        context   = zmq.Context(1)
        publisher = context.socket(zmq.PUB)
        publisher.bind(self.url)

        while True:
            xOpt = self.xOptQueue.get()
            # empty the queue
            while not self.xOptQueue.empty():
                xOpt = self.xOptQueue.get()
                myiter += 1
            traj = trajectory.Trajectory(self.ocp, xOpt)
            for callbackFun, zeromqChannel in self.callbackFunsWithChannels:
                mcStr = callbackFun(traj,myiter,self.ocp)
                publisher.send_multipart([zeromqChannel, mcStr])

def trajectoryCallback(toProto,protoTraj,showAllPoints=False):
    def callback(traj,myiter,ocp):
        protos = []
        for k in range(0,ocp.nk):
            for nicpIdx in range(0,ocp.nicp):
                if showAllPoints:
                    degIdxRange = range(ocp.deg+1)
                else:
                    degIdxRange = [0]
                for degIdx in degIdxRange:
                    lookup = lambda name: traj.lookup(name,timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)
                    protos.append( toProto(lookup) )
        mc = protoTraj()
        mc.traj.extend(list(protos))

        return mc.SerializeToString()
    return callback
