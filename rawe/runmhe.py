import numpy as np

import casadi as C

from collocation import LagrangePoly
import models
from newton.newton import Newton
from newton.nmhe import Nmhe

if __name__ == '__main__':
    print "creating model"
    dae = models.pendulum2()
#    dae.addP('endTime')

    nk = 20
    nicp = 1
    deg = 4
    newton = Newton(LagrangePoly,dae,nk,nicp,deg,'RADAU')
    endTime = 0.025*nk
    newton.setupStuff(endTime)

    r = 0.3
    x0 = [r,0,0,0]

    xTraj = [np.array(x0)]
    uTraj = []
    p = np.array([0.3])
    
    # make a trajectory
    for k in range(nk):
        uTraj.append(np.array([5*C.sin(k/float(nk))]))
#        uTraj.append(np.array([0]))
        newton.isolver.setInput(xTraj[-1],0)
        newton.isolver.setInput(uTraj[-1],1)
        newton.isolver.setInput(p[-1],2)
        newton.isolver.output().set(1)
        newton.isolver.evaluate()
        xTraj.append(np.array(newton.isolver.output(1)).squeeze())

#    for k in range(nk):
#        print xTraj[k],uTraj[k]
#    print xTraj[-1]

    # setup MHE
    nmhe = Nmhe(dae,nk)

    # bounds
    r = 0.3
    nmhe.bound('x',(-2*r,2*r))
    nmhe.bound('z',(-2*r,2*r))
    nmhe.bound('dx',(-50,50))
    nmhe.bound('dz',(-50,50))
    nmhe.bound('m',(0.28,0.32))

    # boundary conditions
#    nmhe.bound('x',(r,r),timestep=0)
#    nmhe.bound('z',(0,0),timestep=0)
#    nmhe.bound('dx',(0,0),timestep=0)
#    nmhe.bound('dz',(0,0),timestep=0)

    # constrain invariants
    nmhe.constrain(nmhe.lookup('c',timestep=0),'==',0)
    nmhe.constrain(nmhe.lookup('cdot',timestep=0),'==',0)

    # initial guess
    nmhe.guess('m',0.3)
    for k in range(nk+1):
        nmhe.guess('x', xTraj[k][0],timestep=k)
        nmhe.guess('z', xTraj[k][1]+0.1,timestep=k)
        nmhe.guess('dx',xTraj[k][2],timestep=k)
        nmhe.guess('dz',xTraj[k][3],timestep=k)

    # make objective
    obj = 0
    for k in range(nk+1):
        nmhe.addGaussNewtonObjF(nmhe('x',timestep=k) - xTraj[k][0])
        nmhe.addGaussNewtonObjF(nmhe('z',timestep=k) - xTraj[k][1])
        nmhe.addGaussNewtonObjF(nmhe('dx',timestep=k) - xTraj[k][2])
        nmhe.addGaussNewtonObjF(nmhe('dz',timestep=k) - xTraj[k][3])

#        obj += (nmhe('x',timestep=k) - xTraj[k][0])**2
#        obj += (nmhe('z',timestep=k) - xTraj[k][1])**2
#        obj += (nmhe('dx',timestep=k) - xTraj[k][2])**2
#        obj += (nmhe('dz',timestep=k) - xTraj[k][3])**2
        
#        obj += 1e-8*nmhe('x',timestep=k)**2
#        obj += 1e-8*nmhe('z',timestep=k)**2
#        obj += 1e-8*nmhe('dx',timestep=k)**2
#        obj += 1e-8*nmhe('dz',timestep=k)**2
#    obj += 1e-8*nmhe('m')**2


#    nmhe.constrain(nmhe('dz',timestep=0)**2,'<=',1000)

    nmhe.setObj(obj)
    uTraj = C.DMatrix(np.concatenate(uTraj))
    nmhe.makeSolver()
    nmhe.runSolver(uTraj)



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
