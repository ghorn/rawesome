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

import rawe
from rawe.joy import Joy
from rawekite.crosswind_steady_state import getSteadyState
import carousel_dae
import autogen

class Communicator(object):
    def __init__(self):
        self.context   = zmq.Context(1)
        self.publisher = self.context.socket(zmq.PUB)
        self.publisher.bind("tcp://*:5563")
#        self.fOutputs = fOutputs
#        self.outputNames = outputNames

    def send_kite(self, x_sim, u_sim, p_sim, outs_sim, other_messages=None):
        '''
        take in x/u/p/outs dictionaries and a list of string messages
        send a packet over zmq
        '''
        traj_pb = autogen.carousel_pb2.Trajectory()
        def lookup(name):
            for some_dict in [x_sim, u_sim, p_sim, outs_sim]:
                try:
                    return some_dict[name]
                except:
                    pass
        dae_pb = autogen.tocarouselProto.toProto(lookup)
        traj_pb.traj.extend([dae_pb])
        if other_messages is not None:
            for other_message in other_messages:
                traj_pb.messages.append(other_message)
        self.publisher.send_multipart(['carousel sim',
                                       dae_pb.SerializeToString()])
        self.publisher.send_multipart(['carousel', traj_pb.SerializeToString()])

    def close(self):
        '''
        close the zeromq publisher and context
        '''
        self.publisher.close()
        self.context.term()

def run_sim():
    # create the model
    dae = carousel_dae.makeDae()

    # compute the steady state
    steady_state, _ = getSteadyState(dae, 100, 20)

    # create the sim
    dt = 0.01
    sim = rawe.sim.Sim(dae, dt)
    communicator = Communicator()
    joy = Joy()

    # set the initial state from steady_state
    sim_x = {}
    for name in dae.xNames():
        sim_x[name] = steady_state[name]
    sim_u = {}
    for name in dae.uNames():
        sim_u[name] = steady_state[name]
    sim_p = {}
    for name in dae.pNames():
        sim_p[name] = steady_state[name]

    # set up the sim timer
    timer = rawe.sim.Timer(dt)
    timer.start()
    # loop through and simulate
    # if there's an error close the communicator and throw exception
    sim_time = 0.0
    try:
        while True:
            # sleep for dt
            timer.sleep()
            print "running: "+str(sim_time)
            # send message to visualizer/plotter
            communicator.send_kite(sim_x, sim_u, sim_p,
                                   sim.getOutputs(sim_x, sim_u, sim_p))
            # try to take a simulation step of dt
            try:
                js = joy.getAll()
                rudder = js['axes'][0]*numpy.radians(5)
                aileron = js['axes'][3]*numpy.radians(5)
                elevator = js['axes'][4]*numpy.radians(5)
                #print "rudder: ",rudder
                #print "aileron: ",aileron
                #print "elevator: ",elevator
                sim_x['rudder'] = rudder
                sim_x['aileron'] = aileron
                sim_x['elevator'] = elevator
                sim_x = sim.step(sim_x, sim_u, sim_p)
            except RuntimeError:
                # problem simulating, close the communicator
                communicator.close()
                raise Exception('OH NOES, IDAS CHOKED')
            sim_time += dt
    except KeyboardInterrupt:
        print "closing..."
        communicator.close()

if __name__ == '__main__':
    run_sim()
