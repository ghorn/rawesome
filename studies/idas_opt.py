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
#import time
#import os

import numpy
from numpy import pi
import copy

import casadi as C

import kite_pb2
import kiteproto

from multipleShooting import MultipleShootingStage
import models

#tc0 = 2*389.970797939731

x0 = C.DMatrix( [ 1.154244772411
                , -0.103540608242
                , -0.347959211327
                , 0.124930983341
                , 0.991534857363
                , 0.035367725910
                , 0.316039689643
                , -0.073559821379
                , 0.945889986864
                , 0.940484536806
                , -0.106993361072
                , -0.322554269411
                , 0.000000000000
                , 0.000000000000
                , 0.000000000000
                , 0.137035790811
                , 3.664945343102
                , -1.249768772258
                , 0.000000000000
                , 3.874600000000
                ])
x0=C.veccat([x0,C.sqrt(C.sumAll(x0[0:2]*x0[0:2])),0])
rArm = 1.085 #(dixit Kurt)
zt = -0.03

def main():
    nSteps = 15

    print "creating model"
    dae = models.carousel(zt,rArm,nSteps)

    print "setting up OCP"
    ocp = MultipleShootingStage(dae, nSteps)

    # make the integrator
    print "setting up dynamics constraints"
    integratorOptions = [ ("reltol",1e-7)
                        , ("abstol",1e-9)
                        , ("t0",0)
                        , ("tf",1)
                        , ('name','integrator')
                        , ("linear_solver_creator",C.CSparse)
#                        , ("linear_solver_creator",C.LapackLUDense)
                        , ("linear_solver","user_defined")
                        ]
    ocp.setIdasIntegrator(integratorOptions)

    # constrain invariants
    def invariantErrs():
        dcm = C.horzcat( [ C.veccat([dae.x('e11'), dae.x('e21'), dae.x('e31')])
                         , C.veccat([dae.x('e12'), dae.x('e22'), dae.x('e32')])
                         , C.veccat([dae.x('e13'), dae.x('e23'), dae.x('e33')])
                         ] ).trans()
        err = C.mul(dcm.trans(), dcm)
        dcmErr = C.veccat([ err[0,0]-1, err[1,1]-1, err[2,2]-1, err[0,1], err[0,2], err[1,2] ])
        f = C.SXFunction( [dae.xVec(),dae.uVec(),dae.pVec()]
                        , [dae.output('c'),dae.output('cdot'),dcmErr]
                        )
        f.setOption('name','invariant errors')
        f.init()
        return f

    [c0,cdot0,dcmError0] = invariantErrs().call([ocp.states[:,0],ocp.actions[:,0],ocp.params])
    ocp.addConstraint(c0,'==',0)
    ocp.addConstraint(cdot0,'==',0)
    ocp.addConstraint(dcmError0,'==',0)

    # bounds
    ocp.setBound('aileron',(-0.04,0.04))
    ocp.setBound('elevator',(-0.1,0.1))

    ocp.setBound('x',(0,4))
    ocp.setBound('y',(-3,3))
    ocp.setBound('z',(-2,3))
    ocp.setBound('r',(1,2))
    ocp.setBound('dr',(-1,1))
    ocp.setBound('ddr',(0,0))
    ocp.setBound('r',(1.2,1.2),timestep=0)

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.setBound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.setBound(d,(-50,50))

    for w in ['w1','w2','w3']:
        ocp.setBound(w,(-8*pi,8*pi))

    ocp.setBound('delta',(-0.01,1.01*2*pi))
    ocp.setBound('ddelta',(-pi/4,8*pi))
    ocp.setBound('tc',(-200,600))
#    ocp.setBound('tc',(389.970797939731,389.970797939731))
    ocp.setBound('endTime',(0.5,2.0))
#    ocp.setBound('endTime',(1.6336935276077966,1.6336935276077966))
    ocp.setBound('w0',(10,10))

    # boundary conditions
    ocp.setBound('delta',(0,0),timestep=0)
    ocp.setBound('delta',(2*pi,2*pi),timestep=nSteps-1)

    # make it periodic
    ocp.addConstraint(ocp.states[:18,0],'==',ocp.states[:18,-1])
    ocp.addConstraint(ocp.states[19,0],'==',ocp.states[19,-1])
    ocp.addConstraint(ocp.actions[:,0],'==',ocp.actions[:,-1])

    # make the solver
    # objective function
    tc0 = 390
    obj = (C.sumAll(ocp.actions[0:2,:]*ocp.actions[0:2,:]) + 1e-10*C.sumAll((ocp.actions[2,:]-tc0)*(ocp.actions[2,:]-tc0)))*ocp.lookup('endTime')
    ocp.setObjective(obj)

    # zero mq setup
    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    # callback function
    class MyCallback:
      def __init__(self):
        self.iter = 0
      def __call__(self,f,*args):
          self.iter = self.iter + 1
          xOpt = f.input(C.NLP_X_OPT)

          xs,us,p = ocp.getTimestepsFromDvs(xOpt)
          kiteProtos = [kiteproto.toKiteProto(xs[k],us[k],p,rArm,zt) for k in range(0,nSteps)]

          mc = kite_pb2.KiteOpt()
          mc.horizon.extend(list(kiteProtos))

          xup = ocp.devectorize(xOpt)
          mc.messages.append("endTime: "+str(xup['endTime']))
          mc.messages.append("w0: "+str(xup['w0']))
          mc.messages.append("iter: "+str(self.iter))

          publisher.send_multipart(["multi-carousel", mc.SerializeToString()])

    def makeCallback():
        nd = ocp.getDesignVars().size()
        nc = ocp._constraints.getG().size()
        c = C.PyFunction( MyCallback(), C.nlpsolverOut(x_opt=C.sp_dense(nd,1), cost=C.sp_dense(1,1), lambda_x=C.sp_dense(nd,1), lambda_g = C.sp_dense(nc,1), g = C.sp_dense(nc,1) ), [C.sp_dense(1,1)] )
        c.init()
        return c


    # solver
    solverOptions = [ ("iteration_callback",makeCallback())
                    , ("linear_solver","ma57")
#                    , ("max_iter",5)
                    ]

    ocp.setSolver( C.IpoptSolver, solverOptions=solverOptions )
    #ocp.setBounds()

    # initial conditions
    ocp.setXGuess(x0)
    for k in range(0,nSteps):
        val = 2*pi*k/(nSteps-1)
        ocp.setGuess('delta',val,timestep=k,quiet=True)

    ocp.setGuess('aileron',0)
    ocp.setGuess('elevator',0)
    ocp.setGuess('tc',389.970797939731)
    ocp.setGuess('endTime',1.6336935276077966)

    ocp.setGuess('ddr',0)
    ocp.setGuess('w0',0)

    ocp.solve()

if __name__=='__main__':
    main()
