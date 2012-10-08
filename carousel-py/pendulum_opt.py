import zmq

import numpy
from numpy import pi
import copy

import kite_pb2
import casadi as C

from ocputils import MultipleShootingInterval
import pendulum_model

def main():
    nSteps = 15

    print "creating model"
    dae = pendulum_model.pendulum_model(nSteps=nSteps)

    print "setting up OCP"
    ocp = MultipleShootingInterval(dae, nSteps)
    
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
          
          po = kite_pb2.PendulumOpt()
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
    ocp.solve()

if __name__=='__main__':
    main()
