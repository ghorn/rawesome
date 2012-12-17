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
    endTime = 0.5
    ocp.setupCollocation(endTime)

    ffcn = ocp._makeResidualFun()
    x0 = C.ssym('x',dae.xVec().size())
    X = C.ssym('x',dae.xVec().size(),ocp.deg)
    Z = C.ssym('z',dae.zVec().size(),ocp.deg)
    U = C.ssym('u',dae.uVec().size())
    P = C.ssym('p',dae.pVec().size())
    constraints = []
    ############################################################
    ndiff = dae.xVec().size()
    for j in range(1,ocp.deg+1):
        # Get an expression for the state derivative at the collocation point
        xp_jk = ocp.lagrangePoly.lDotAtTauRoot[j,0]*x0
        for j2 in range (1,ocp.deg+1):
            # get the time derivative of the differential states (eq 10.19b)
            xp_jk += ocp.lagrangePoly.lDotAtTauRoot[j,j2]*X[:,j2-1] #ocp.xVec(k,nicpIdx=i,degIdx=j2)
        # Add collocation equations to the NLP
        [fk] = ffcn.eval([xp_jk/ocp.h,
                          X[:,j-1],
                          Z[:,j-1],
                          U,
                          P])
        
        # impose system dynamics (for the differential states (eq 10.19b))
        constraints.append(fk[:ndiff])
#        ocp.constrain(fk[:ndiff],'==',0,tag=
#                       "system dynamics, differential states, kIdx: %d,nicpIdx: %d, degIdx: %d" % (k,i,j))
        
        # impose system dynamics (for the algebraic states (eq 10.19b))
        constraints.append(fk[ndiff:])
#        ocp.constrain(fk[ndiff:],'==',0,tag=
#                       "system dynamics, algebraic states, kIdx: %d,nicpIdx: %d, degIdx: %d" % (k,i,j))
        
    # Get an expression for the state at the end of the finite element
    xf = ocp.lagrangePoly.lAtOne[0]*x0
    for j in range(1,ocp.deg+1):
        xf += ocp.lagrangePoly.lAtOne[j]*X[:,j-1]

#    ifcnSx = C.SXFunction([C.veccat([xUnknown,zUnknown,xf]),x0,u,m],[constraints,xf])
    ifcnSx = C.SXFunction([C.veccat([X,Z]),x0],[C.veccat(constraints),xf])
    ifcnSx.init()

    f = C.SXFunction([C.veccat([X,Z]),x0],[0])
    f.init()
    g = C.SXFunction([C.veccat([X,Z]),x0],[C.veccat(constraints)])
    g.init()

    solver = C.IpoptSolver(f,g)
    solver.setOption("parametric",True)
    solver.init()

    solver.input(C.NLP_LBG).set(0)
    solver.input(C.NLP_UBG).set(0)
#    solver.setInput(0,C.NLP_LBG)
#    solver.setInput(0,C.NLP_UBG)

    r = 0.3
    x0_ = C.veccat([r,0,0,0])
    z0_ = C.veccat([0])
    
    XInit = C.horzcat([x0_]*ocp.deg)
    ZInit = C.horzcat([z0_]*ocp.deg)
    
    solver.setInput([0.3,0,0,0],C.NLP_P)
    solver.setInput(C.veccat([XInit,ZInit]), C.NLP_X_INIT)
    solver.evaluate()
    xopt = solver.output(C.NLP_X_OPT)
    x = xopt[:X.size()]
    z = xopt[X.size():]

    x = x.reshape((4,ocp.deg))
    print x

#    print C.veccat(constraints),xf
    import sys;sys.exit()
#    ifcn = C.KinsolSolver(ifcnSx)
    ifcn = C.NLPImplicitSolver(ifcnSx)
    ifcn.setOption("nlp_solver",C.IpoptSolver)
    ifcn.init()

    r = 0.3
    x0 = [r,0,0,0]
#    ifcn.input(0).set(x0)
    ifcn.setInput(x0,0)
#    ifcn.setInput(0,1) # torque
#    ifcn.setInput(0.3,2) # mass
    ifcn.evaluate()
    
#    C.SXFunction([xUnknown,
#    ifun = C.MXFunction( [V[4:],V[0:4]], [constraints,ocp.xVec(1)] )

#    sx_unknowns = ssym("unknowns",unknowns.shape)
#    sx_knowns = [ssym("knowns",k.shape) for k in knowns]
#
#    [5,6,7,8,10,11,12,13,15,16,17,18,20,21,22,23,9,14,19,24,]
#
#    fun.eval([sx_unknowns+sx_knowns])
#
#    fun.init()

#    ifun = C.MXFunction( [V[4:],V[0:4]], [constraints,ocp.xVec(1)] )
    ifcn.jacobian(1,0)
