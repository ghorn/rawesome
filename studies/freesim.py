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
import time
import os
import simutils

import casadi as C

import kiteproto
import models
import joy

#tc0 = 2*389.970797939731

x0 = C.DMatrix( [ -10 # position
                , 0
                , 15
                , 1   # DCM
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 1
                , 20.0 # vel
                , 0.0
                , 0.0
                , 0   # ang vel
                , 0
                , 0
                ])

ts = 0.02
sim = simutils.Sim(ts=ts, sloMoFactor=4, state0=simutils.SimState(pdOn = False, x = x0))

if __name__=='__main__':
    print "creating model"
    (ode, others, outputs) = models.free()
    ode.init()

    print "creating outputs function"
    outputNames = outputs.keys()
    fOutputs = C.SXFunction([others['xVec'],C.veccat([others['uVec'],others['pVec']])],[outputs[n] for n in outputNames])
    fOutputs.setOption('name','fOutputs')
    fOutputs.init()

    print "creating communicator"
    communicator = simutils.Communicator(fOutputs,outputNames)

    print "creating integrator"
    f = C.CVodesIntegrator(ode)
    f.setOption("reltol",1e-5)
    f.setOption("abstol",1e-7)
    f.setOption("t0",0)
    f.setOption("tf",0.02)
    f.setOption('name','integrator')
#    f.setOption("linear_solver_creator",C.CSparse)
#    f.setOption("linear_solver","user_defined")
#    f.setOption("monitor",["res"])
    f.init()

    def advanceState():
        js = sim.handleInput()

        # saving/loading
        fstButton=13
        for k in range(0,4):
            if k+fstButton+4 in js['buttonsDown'] and k in sim._saves:
                print "loading save #"+str(k+1)
                sim.load(k)
                sim.default = k
            if k+fstButton in js['buttonsDown']:
                if k==0:
                    print "can't override save0"
                else:
                    print "making save #"+str(k+1)
                    sim.save(k)
                    sim.default=k
        if 5 in js['buttonsDown']:
            k = sim.default
            print "loading save #"+str(k+1)
            sim.load(k)

        # play replay
        if 3 in js['buttonsDown']:
            sim.playReplay(communicator)
            sim.loadDefault()

        aileron = -js['axes'][0]*0.05
        elevator =  js['axes'][1]*0.2
        rudder = -js['axes'][2]*0.15
        tc = 600*(1 - js['axes'][6])
        wind_x = 5*(1-js['axes'][7])
        wind_x = 0

        x = sim.getCurrentState().x
        u = C.DMatrix([tc,aileron,elevator,rudder])
        p = C.DMatrix([wind_x])

        f.setInput(x,C.INTEGRATOR_X0)
        f.setInput(C.veccat([u,p]),C.INTEGRATOR_P)
        f.evaluate()

        xNext = C.DMatrix(f.output())
        return ((x, u, p), xNext)

    x = x0
    print "simulating..."
    try:
        while True:
            t0 = time.time()
            try:
                ((x,u,p), xNext) = advanceState()
                sim.currentState.log(x,u,p)
                sim.currentState.x = xNext
                if len(sim._saves[sim.default]._log)==0:
                    sim.save(sim.default)
            except RuntimeError:
                sim.loadDefault()
                x,u,p = sim.getCurrentState()._log[-1]
                pass
            communicator.sendKite(sim,(x,u,p))

            deltaTime = (t0 + sim.tsSimStep*sim.sloMoFactor) - time.time()
            if deltaTime > 0:
                time.sleep(deltaTime)
    except KeyboardInterrupt:
        print "closing..."
        publisher.close()
        context.term()
        pass
