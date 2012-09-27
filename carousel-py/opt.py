#import time
#import os

import numpy
from numpy import pi
import copy

import casadi as C

import ocp 
import model

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

r = 1.2

def toProto(x,u):
    cs = kite_pb2.CarouselState()

    cs.kiteXyz.x = x.at(0)
    cs.kiteXyz.y = x.at(1)
    cs.kiteXyz.z = x.at(2)

    cs.kiteDcm.r11 = x.at(3)
    cs.kiteDcm.r12 = x.at(4)
    cs.kiteDcm.r13 = x.at(5)

    cs.kiteDcm.r21 = x.at(6)
    cs.kiteDcm.r22 = x.at(7)
    cs.kiteDcm.r23 = x.at(8)

    cs.kiteDcm.r31 = x.at(9)
    cs.kiteDcm.r32 = x.at(10)
    cs.kiteDcm.r33 = x.at(11)

    cs.delta = x.at(18)
    cs.ddelta = x.at(19)
    
    cs.tc = u.at(0)
    cs.u0 = u.at(1)
    cs.u1 = u.at(2)
    return cs

        

def main():
    nSteps = 10
    endTime = C.ssym('endTime')

    print "creating model"
    (dae, others) = model.model((endTime,nSteps))
    dae.init()

#    # compile model code
#    print "generating model code"
#    t0 = time.time()
#    dae.generateCode("dae.c")
#    print "took "+str(time.time()-t0)+" seconds to generate code"
#    t0 = time.time()
#    os.system("gcc -fPIC -O2 -shared dae.c -o dae.so")
#    print "took "+str(time.time()-t0)+" seconds to compile code"
#    dae_ext = C.ExternalFunction("./dae.so")
#    dae_ext.init()
#    dae = dae_ext

    nStates = others['xVec'].size()
    nActions = others['uVec'].size()
    nParams = others['pVec'].size()

    assert(nStates==dae.inputSX(C.DAE_X).size())
    assert(nActions+nParams==dae.inputSX(C.DAE_P).size())
    
    # make the integrator
    print "creating integrator"
    integrator = C.IdasIntegrator(dae)
    integrator.setOption("reltol",1e-6)
    integrator.setOption("abstol",1e-8)
    integrator.setOption("t0",0)
    integrator.setOption("tf",1)
    integrator.setOption('name','integrator')
    integrator.init()

    # make the OCP
    print "setting up OCP"
    states = C.msym("x" ,nStates,nSteps)
    actions = C.msym("u",nActions,nSteps)
    params = C.msym("p",nParams)
    constraints = ocp.Constraints()

    constraints.addDynamicsConstraints(integrator, states, actions, params)

    # constrain invariants
    def invariantErrs():
        dcm = C.horzcat( [ C.veccat([others['xDict']['e11'], others['xDict']['e21'], others['xDict']['e31']])
                         , C.veccat([others['xDict']['e12'], others['xDict']['e22'], others['xDict']['e32']])
                         , C.veccat([others['xDict']['e13'], others['xDict']['e23'], others['xDict']['e33']])
                         ] ).trans()
        err = C.mul(dcm.trans(), dcm)
        dcmErr = C.veccat([ err[0,0]-1, err[1,1]-1, err[2,2]-1, err[0,1], err[0,2], err[1,2] ])
        f = C.SXFunction( [others['xVec'],others['uVec'],others['pVec']]
                        , [others['c'],others['cdot'],dcmErr]
                        )
        f.setOption('name','invariant errors')
        f.init()
        return f
    
    f = invariantErrs()
    [c0,cdot0,dcmError0] = f.call([states[:,0],actions[:,0],params])
    constraints.add(c0,'==',0)
    constraints.add(cdot0,'==',0)
    constraints.add(dcmError0,'==',0)

    # make it periodic
    constraints.add(states[:,0],'==',states[:,-1])
    constraints.add(actions[:,0],'==',actions[:,-1])

    # bounds
    bounds = ocp.Bounds(others['xNames'], others['uNames'], others['pNames'], nSteps)
    bounds.setBound('u1',(-0.04,0.04))
    bounds.setBound('u2',(-0.1,0.1))
    
    bounds.setBound('x',(r-2,r+2))
    bounds.setBound('y',(-3,3))
    bounds.setBound('z',(-3,3))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        bounds.setBound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        bounds.setBound(d,(-50,50))

    for w in ['w1','w2','w3']:
        bounds.setBound(w,(-4*pi,4*pi))

    bounds.setBound('delta',(-0.01,1.01*2*pi))
    bounds.setBound('ddelta',(-pi/4,8*pi))
    bounds.setBound('tc',(-200,600))

    bounds.setBound('endTime',(0.5,5))
    
    bounds.setBound('wind_x',(0,0))

    # boundary conditions
    bounds.setBound('delta',(0,0),0)
    bounds.setBound('delta',(2*pi,2*pi),nSteps-1)

    # make the solver
    designVars = C.veccat( [C.flatten(states), C.flatten(actions), C.flatten(params)] )
    
    # objective function
    obj = C.sumAll(actions*actions)
    f = C.MXFunction([designVars], [obj])

    # constraint function
    g = C.MXFunction([designVars], [constraints.getG()])

    # solver
    solver = C.IpoptSolver(f, g)
    solver.init()

    solver.setInput(constraints.getLb(), C.NLP_LBG)
    solver.setInput(constraints.getUb(), C.NLP_UBG)

    lb,ub = bounds.get()
    solver.setInput(lb, C.NLP_LBX)
    solver.setInput(ub, C.NLP_UBX)

    # initial conditions
    guess = ocp.InitialGuess(others['xNames'], others['uNames'], others['pNames'], nSteps)
    guess.setXVec(x0)
    for k in range(0,nSteps):
        val = 2*pi*k/(nSteps-1)
        guess.setGuess('delta',val,timestep=k,quiet=True)

    guess.setGuess('u1',0)
    guess.setGuess('u2',0)
    guess.setGuess('tc',389.970797939731)
    guess.setGuess('endTime',1.6336935276077966)
    
    guess.setGuess('wind_x',0)
    
    solver.setInput(guess.get(), C.NLP_X_INIT)

    solver.solve()

if __name__=='__main__':
    main()
