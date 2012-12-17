import casadi as C

import models
from collocation import Coll

if __name__ == '__main__':
    print "creating model"
    dae = models.pendulum2()
#    dae.addP('endTime')

#    x = dae.xVec()
#    z = dae.zVec()
#    u = dae.uVec()
#    p = dae.pVec()
#    xDot = C.veccat([dae.ddt(name) for name in dae.xNames()])

    ocp = Coll(dae, nk=1, nicp=1, collPoly="RADAU", deg=4)
    endTime = 0.01
    ocp.setupCollocation(endTime)

    ############################################################

    V = ocp._dvMap.vectorize()
    constraints = ocp._constraints.getG()
    
    m = C.MXFunction([V],[constraints,
                          C.veccat([ocp.xVec(0,degIdx=k) for k in range(1,ocp.deg+1)]),
                          C.veccat([ocp.zVec(0,degIdx=k) for k in range(1,ocp.deg+1)]),
                          ocp.xVec(0,degIdx=0),
                          ocp.xVec(1,degIdx=0),
                          ocp('torque',timestep=0),
                          ocp('m')])
    m.init()

    s = C.SXFunction(m)
    s.init()
    sxV = C.ssym("sxV",V.shape)
#    [constraints,xUnknown,zUnknown,x0,xf] = s.eval([sxV])
    [constraints,xUnknown,zUnknown,x0,xf,u,m] = s.eval([sxV])
#    print xUnknown
#    print zUnknown
#    print x0
#    print xf

    ifcnSx = C.SXFunction([C.veccat([xUnknown,zUnknown,xf]),x0,u,m],[constraints,xf])
#    ifcnSx = C.SXFunction([C.veccat([xUnknown,zUnknown,xf]),x0],[constraints,xf])
    ifcnSx.init()
    
#    ifcn = C.KinsolSolver(ifcnSx)
    ifcn = C.NLPImplicitSolver(ifcnSx)
    ifcn.setOption("nlp_solver",C.IpoptSolver)
    ifcn.setOption("linear_solver",C.CSparse)
    ifcn.init()

    r = 0.3
    x0 = [r,0,0,0]
    ifcn.input(0).set(x0)
#    ifcn.setInput(x0,0)
    ifcn.setInput(0,1) # torque
    ifcn.setInput(0.3,2) # mass
    ifcn.evaluate()

    j = ifcn.jacobian(0,1)
    j.init()

    j.setInput(x0)
    j.evaluate()
    print j.output(0)
    print j.output(1)
