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

import matplotlib.pyplot as plt
import zmq
import numpy
from numpy import pi
import copy

import casadi as C

import rawe

from multipleShooting import MultipleShootingStage
import models

def main():
    nSteps = 15

    print "creating model"
    dae = models.pendulum(nSteps=nSteps)

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
        f = C.SXFunction( [dae.xVec(),dae.uVec(),dae.pVec()]
                        , [dae.output('invariants')]
                        )
        f.setOption('name','invariant errors')
        f.init()
        return f

    [c0] = invariantErrs().call([ocp.states[:,0],ocp.actions[:,0],ocp.params])
    ocp.addConstraint(c0,'==',0)

    # bounds
    r = 0.3
    ocp.setBound('x',(-0.5,0.5))
    ocp.setBound('z',(-0.5,0.5))
    ocp.setBound('dx',(-5,5))
    ocp.setBound('dz',(-5,5))
    ocp.setBound('torque',(-50,50))
    ocp.setBound('m',(0.3,0.5))
    ocp.setBound('endTime',(0.01,5.0))

    # boundary conditions
    ocp.setBound('x',(r,r),timestep=0)
    ocp.setBound('z',(0,0),timestep=0)
    ocp.setBound('x',(0,0),timestep=nSteps-1)
    ocp.setBound('z',(-r*1.5,-r/2),timestep=nSteps-1)
    ocp.setBound('dx',(0,0),timestep=0)
    ocp.setBound('dz',(0,0),timestep=0)
    ocp.setBound('dx',(0,0),timestep=nSteps-1)
    ocp.setBound('dz',(0,0),timestep=nSteps-1)

    # make the solver
    ocp.setObjective(ocp.lookup('endTime'))

    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    # callback function
    class MyCallback:
      def __init__(self):
        self.iter = 0
      def __call__(self,f,*args):
          xOpt = f.input(C.NLP_X_OPT)
          self.iter = self.iter + 1
          xup = ocp.devectorize(xOpt)

          po = rawe.kite_pb2.PendulumOpt()
          po.x.extend(list(xup['x']))
          po.z.extend(list(xup['z']))
          po.endTime = xup['endTime']
          po.iters = self.iter
          publisher.send_multipart(["pendulum-opt", po.SerializeToString()])

    def makeCallback():
        nd = ocp.getDesignVars().size()
        nc = ocp._constraints.getG().size()
        c = C.PyFunction( MyCallback(), C.nlpsolverOut(x_opt=C.sp_dense(nd,1), cost=C.sp_dense(1,1), lambda_x=C.sp_dense(nd,1), lambda_g = C.sp_dense(nc,1), g = C.sp_dense(nc,1) ), [C.sp_dense(1,1)] )
        c.init()
        return c

    # solver
    solverOptions = [ ("iteration_callback",makeCallback())
#                    , ("derivative_test","first-order")
                    , ("linear_solver","ma57")
                    ]
    constraintFunOptions = [('numeric_jacobian',False)]
    ocp.setSolver( C.IpoptSolver, solverOptions=solverOptions, constraintFunOptions=constraintFunOptions )

    # initial conditions
    ocp.setXGuess([r,0,0,0])
    ocp.setGuess('torque',0)
    ocp.setGuess('m',0.4)
    ocp.setGuess('endTime',0.5)

    #ocp.setBounds()

    # solve!
    opt = ocp.devectorize(ocp.solve())

    # Plot the results
    plt.figure(1)
    plt.clf()
    time = numpy.linspace(0,opt['endTime'],opt['x'].size())
    plt.plot(time, opt['x'], '--')
    plt.plot(time, opt['z'], '-')
    plt.plot(time, opt['dx'], '-.')
    plt.plot(time, opt['dz'], '--')

    time = numpy.linspace(0,opt['endTime'],opt['torque'].size())
    plt.plot(time,opt['torque']/20,'--')

    plt.title("pendulum swingup optimization")
    plt.xlabel('time')
    plt.legend(['x','z','vx','vz','torque'])
    plt.grid()
    plt.show()


if __name__=='__main__':
    main()
