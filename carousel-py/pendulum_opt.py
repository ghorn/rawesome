#import time
#import os

import numpy
from numpy import pi
import copy

import casadi as C

import ocp 
import pendulum_model

def main():
    nSteps = 15
    endTime = C.ssym('endTime')

    print "creating model"
    (dae, others) = pendulum_model.pendulum_model((endTime,nSteps))
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
        f = C.SXFunction( [others['xVec'],others['uVec'],others['pVec']]
                        , [C.veccat(others['c'])]
                        )
        f.setOption('name','invariant errors')
        f.init()
        return f
    
    [c0] = invariantErrs().call([states[:,0],actions[:,0],params])
    constraints.add(c0,'==',0)

    # bounds
    bounds = ocp.Bounds(others['xNames'], others['uNames'], others['pNames'], nSteps)
    r = 0.3
    bounds.setBound('x',(-0.5,0.5))
    bounds.setBound('z',(-0.5,0.5))
    bounds.setBound('dx',(-5,5))
    bounds.setBound('dz',(-5,5))
    bounds.setBound('torque',(-15,15))
    bounds.setBound('m',(0.3,0.5))
    bounds.setBound('endTime',(0.2,2.0))

    # boundary conditions
    bounds.setBound('x',(r,r),timestep=0)
    bounds.setBound('z',(0,0),timestep=0)
    bounds.setBound('x',(-0.01,0.01),timestep=nSteps-1)
    bounds.setBound('z',(-r*1.01,-r*0.99),timestep=nSteps-1)
#    bounds.setBound('x',(0,0),timestep=nSteps-1)
#    bounds.setBound('z',(-r,-r),timestep=nSteps-1)
    bounds.setBound('dx',(0,0),timestep=0)
    bounds.setBound('dz',(0,0),timestep=0)
    bounds.setBound('dx',(0,0),timestep=nSteps-1)
    bounds.setBound('dz',(0,0),timestep=nSteps-1)

    # make the solver
    designVars = C.veccat( [C.flatten(states), C.flatten(actions), C.flatten(params)] )
    
    # objective function
    obj = C.sumAll(actions*actions)
    f = C.MXFunction([designVars], [obj])

    # constraint function
    g = C.MXFunction([designVars], [constraints.getG()])

    # callback function
    class MyCallback:
      def __init__(self):
        self.iter = 0 
#      def __call__(self,f,*args):
#        print "====Hey, I'm an iteration===="
#        print "X_OPT = ", f.input(C.NLP_X_OPT)
#        print f.getStats()
#        self.iter = self.iter + 1
#        if self.iter > 5:
#          print "5 Iterations, that is quite enough!"
#          f.output(0).set(1) # With this statement you can halt the iterations
      def __call__(self,f,*args):
        x = f.input(C.NLP_X_OPT)
        self.iter = self.iter + 1
        x = bounds.devectorize(x)
        x['x']
        
    def makeCallback():
        nd = designVars.size()
        nc = constraints.getG().size()
        c = C.PyFunction( MyCallback(), C.nlpsolverOut(x_opt=C.sp_dense(nd,1), cost=C.sp_dense(1,1), lambda_x=C.sp_dense(nd,1), lambda_g = C.sp_dense(nc,1), g = C.sp_dense(nc,1) ), [C.sp_dense(1,1)] )
        c.init()
        return c
        

    # solver
    solver = C.IpoptSolver(f, g)
    solver.setOption("iteration_callback",makeCallback())
#    solver.setOption("derivative_test","first-order")
    solver.setOption("linear_solver","ma57")
    solver.init()

    solver.setInput(constraints.getLb(), C.NLP_LBG)
    solver.setInput(constraints.getUb(), C.NLP_UBG)

    lb,ub = bounds.get()
    solver.setInput(lb, C.NLP_LBX)
    solver.setInput(ub, C.NLP_UBX)

    # initial conditions
    guess = ocp.InitialGuess(others['xNames'], others['uNames'], others['pNames'], nSteps)
    guess.setXVec([r,0,0,0])
    guess.setGuess('torque',0)
    guess.setGuess('m',0.4)
    
    guess.setGuess('endTime',0.5)
    
    solver.setInput(guess.get(), C.NLP_X_INIT)

    solver.solve()
    xopt = solver.output(C.NLP_X_OPT)

    f.setInput(xopt,0)
    f.evaluate()
    print f.output(0)

    g.setInput(xopt,0)
    g.evaluate()
    print g.output(0)

if __name__=='__main__':
    main()
