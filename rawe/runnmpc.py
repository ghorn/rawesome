import numpy as np

import casadi as C

from collocation import Coll,LagrangePoly
import models
from newton.newton import Newton
from newton.multipleShooting import Nmpc

if __name__ == '__main__':
    print "creating model"
    dae = models.pendulum2()
    dae.addP('endTime')

    nk = 30
    nicp = 1
    deg = 4

    nmpc = Nmpc(dae,nk)
    
    # constrain invariants
    nmpc.constrain(nmpc.lookup('c',timestep=0),'==',0)
    nmpc.constrain(nmpc.lookup('cdot',timestep=0),'==',0)

    # bounds
    r = 0.3
    nmpc.bound('x',(-2*r,2*r))
    nmpc.bound('z',(-2*r,0.01*r))
    nmpc.bound('dx',(-5,5))
    nmpc.bound('dz',(-5,5))
    nmpc.bound('torque',(-500,500))
    nmpc.bound('m',(0.3,0.3))
    nmpc.bound('endTime',(0.2,3.5))

    # boundary conditions
    nmpc.bound('x',(r,r),timestep=0)
    nmpc.bound('z',(0,0),timestep=0)
    nmpc.bound('x',(0,0),timestep=-1)
    nmpc.bound('z',(-10*r,0.01*r),timestep=-1)
    nmpc.bound('dx',(0,0),timestep=0)
    nmpc.bound('dz',(0,0),timestep=0)
    nmpc.bound('dx',(0,0),timestep=-1)
    nmpc.bound('dz',(-0.5,0.5),timestep=-1)

    # gauss-newton objective
    for k in range(nk):
        nmpc.addGaussNewtonObjF(nmpc('torque',timestep=k))

    # arbitrary objective
    nmpc.setObj(nmpc.lookup('endTime'))

    x0 = [r,0,0,0]

    nmpc.makeSolver()
