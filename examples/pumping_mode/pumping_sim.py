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
import numpy
import casadi as C

import time

import rawe
import rawekite
from rawe.joy import Joy
from rawekite.crosswind_steady_state import getSteadyState
import pumping_dae
import autogen

class Communicator(object):
    def __init__(self):
        self.context   = zmq.Context(1)
        self.publisher = self.context.socket(zmq.PUB)
        self.publisher.bind("tcp://*:5563")
#        self.fOutputs = fOutputs
#        self.outputNames = outputNames

    def sendKite(self,x,u,p,outs,otherMessages=[]):
        pb = autogen.pumping_pb2.Trajectory()
        def lookup(name):
            for d in [x,u,p,outs]:
                try:
                    return d[name]
                except:
                    pass
        dae = autogen.topumpingProto.toProto(lookup)
        pb.traj.extend([dae])
        if len(otherMessages)>0:
            #pb.messages.append("-------------------------")
            for om in otherMessages:
                pb.messages.append(om)
        self.publisher.send_multipart(['pumping sim', dae.SerializeToString()])
#        self.publisher.send_multipart(['pumping trajectory', pb.SerializeToString()])
        self.publisher.send_multipart(['pumping', pb.SerializeToString()])
        
    def close(self):
        self.publisher.close()
        self.context.term()

if __name__=='__main__':
    # create the model
    dae = pumping_dae.makeDae()
    
    # compute the steady state
    steadyState, ssDot = getSteadyState(dae,100,20)
#    print "steady state:"
#    for name in dae.xNames():
#        print name+': '+str(steadyState[name])
#    print "steady state ddt:"
#    print ssDot
#    import sys; sys.exit()
    
    # create the sim
    dt = 0.01
    sim = rawe.sim.Sim(dae,dt)
    communicator = Communicator()
    joy = Joy()

    # set the initial state from steadyState
    x = {}
    for name in dae.xNames():
        x[name] = steadyState[name]
    u = {}
    for name in dae.uNames():
        u[name] = steadyState[name]
    p = {}
    for name in dae.pNames():
        p[name] = steadyState[name]

    # set up the sim timer
    timer = rawe.sim.Timer(dt)
    timer.start()
    # loop through and simulate, if there's an error close the communicator and throw exception
    sim_time = 0.0
    try:
        while True:
            # sleep for dt
            timer.sleep()
            print "running: "+str(sim_time)
            # send message to visualizer/plotter
            outs = sim.getOutputs(x,u,p)
            communicator.sendKite(x,u,p,outs)
            # try to take a simulation step of dt
            try:
                js = joy.getAll()
                rudder = js['axes'][0]*numpy.radians(10)
                aileron = js['axes'][3]*numpy.radians(10)
                elevator = js['axes'][4]*numpy.radians(10)
                #print "rudder: ",rudder
                #print "aileron: ",aileron
                #print "elevator: ",elevator
                x['rudder'] = rudder
                x['aileron'] = aileron
                x['elevator'] = elevator
                x = sim.step(x,u,p)
            except RuntimeError:
                # problem simulating, close the communicator
                communicator.close()
                raise Exception('OH NOES, IDAS CHOKED')
            sim_time += dt
    except KeyboardInterrupt:
        print "closing..."
        communicator.close()
        pass
