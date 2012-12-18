import numpy as np

import casadi as C

from collocation import Coll,LagrangePoly
import models
from newton.newton import Newton,Nmpc

if __name__ == '__main__':
    print "creating model"
    dae = models.pendulum2()
#    dae.addP('endTime')

    nk = 30
    nicp = 1
    deg = 4
    newton = Newton(LagrangePoly,dae,nk,nicp,deg,'RADAU')
    endTime = 0.05*nk
    newton.setupStuff(endTime)

    r = 0.3
    x0 = [r,0,0,0]

    xTraj = [np.array(x0)]
    uTraj = []
    p = np.array([0.3])
    
    # make a trajectory
    for k in range(nk):
        uTraj.append(np.array([5*C.sin(k/float(nk))]))
        newton.isolver.setInput(xTraj[-1],0)
        newton.isolver.setInput(uTraj[-1],1)
        newton.isolver.setInput(p[-1],2)
        if k ==0:
          # NewtonImplicitSolver needs good initialisation
          ifcn = newton.isolver.getF();
          implicitSolver = C.NLPImplicitSolver(ifcn)
          implicitSolver.setOption("nlp_solver",C.IpoptSolver)
          implicitSolver.setOption("nlp_solver_options",{'print_level': 0})
          implicitSolver.setOption("linear_solver",C.CSparse)
          implicitSolver.init()
          for i in range(3):
            implicitSolver.input(i).set(newton.isolver.input(i))
          implicitSolver.evaluate()
          newton.isolver.output().set(implicitSolver.output())
        newton.isolver.evaluate()
        xTraj.append(np.array(newton.isolver.output(1)).squeeze())

    for k in range(nk):
        print xTraj[k],uTraj[k]
    print xTraj[-1]

    # setup MHE
    ocp = Nmpc(dae,nk)

    # make objective
    obj = 0
    for k in range(nk):
        obj += (ocp('x',timestep=k) - xTraj[k][0])**2
        obj += (ocp('z',timestep=k) - xTraj[k][1])**2
        obj += (ocp('dx',timestep=k) - xTraj[k][2])**2
        obj += (ocp('dz',timestep=k) - xTraj[k][3])**2
    print obj



##    newton.isolver.input(0).set(x0)
#    newton.isolver.setInput(x0,0)
#    newton.isolver.setInput(0,1) # torque
#    newton.isolver.setInput(0.3,2) # mass
#    newton.isolver.evaluate()
#    newton.isolver.evaluate()
#
#    j = newton.isolver.jacobian(0,1)
#    j.init()
#    j.setInput(x0,0)
#    j.setInput(0,1)
#    j.setInput(0.3,2)
#    j.evaluate()
#    print j.output(0)
#    print j.output(1)
