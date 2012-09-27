import zmq

import numpy
from numpy import pi
import copy

import kite_pb2
import casadi as C

import ocp 
import pendulum_model

def main():
    nSteps = 15
    endTime = C.ssym('endTime')

    print "creating model"
    (dae, others) = pendulum_model.pendulum_model((endTime,nSteps))
    dae.init()

    nStates = others['xVec'].size()
    nActions = others['uVec'].size()
    nParams = others['pVec'].size()

    assert(nStates==dae.inputSX(C.DAE_X).size())
    assert(nActions+nParams==dae.inputSX(C.DAE_P).size())
    
    # make the integrator
    print "creating integrator"
    integrator = C.IdasIntegrator(dae)
    integrator.setOption("reltol",1e-8)
    integrator.setOption("abstol",1e-10)
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
    bounds.setBound('torque',(-50,50))
    bounds.setBound('m',(0.3,0.5))
    bounds.setBound('endTime',(0.01,5.0))

    # boundary conditions
    bounds.setBound('x',(r,r),timestep=0)
    bounds.setBound('z',(0,0),timestep=0)
    bounds.setBound('x',(0,0),timestep=nSteps-1)
    bounds.setBound('z',(-r*1.5,-r/2),timestep=nSteps-1)
    bounds.setBound('dx',(0,0),timestep=0)
    bounds.setBound('dz',(0,0),timestep=0)
    bounds.setBound('dx',(0,0),timestep=nSteps-1)
    bounds.setBound('dz',(0,0),timestep=nSteps-1)

    # make the solver
    designVars = C.veccat( [C.flatten(states), C.flatten(actions), C.flatten(params)] )
    dvs = ocp.DesignVars((others['xNames'],states), (others['uNames'],actions), (others['pNames'],params), nSteps)
    
    # objective function
    obj = dvs.lookup('endTime')
    f = C.MXFunction([designVars], [obj])

    # constraint function
    g = C.MXFunction([designVars], [constraints.getG()])

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
          xup = bounds.devectorize(xOpt)
          po = kite_pb2.PendulumOpt()
          po.x.extend(list(xup['x']))
          po.z.extend(list(xup['z']))
          po.endTime = xup['endTime']
          po.iters = self.iter
          publisher.send_multipart(["pendulum-opt", po.SerializeToString()])
        
        
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
    
    solver.setInput(guess.vectorize(), C.NLP_X_INIT)

    solver.solve()

if __name__=='__main__':
    main()
