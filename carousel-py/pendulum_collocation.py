import matplotlib.pyplot as plt
import zmq

import numpy
from numpy import pi
import copy

import kite_pb2
import casadi as C

from ocputils import MultipleShootingStage
import pendulum_model
from collocation import Coll

def main():
    nk = 15

    print "creating model"
    dae = pendulum_model.pendulum_model()
    dae.addP('endTime')

    print "setting up OCP"
    ocp = Coll(dae, nk=nk,nicp=1,deg=4)
    
    # constrain invariants
    def invariantErrs():
        f = C.SXFunction( [dae.xVec(),dae.uVec(),dae.pVec()]
                        , [dae.output('invariants')]
                        )
        f.setOption('name','invariant errors')
        f.init()
        return f
    
    [c0] = invariantErrs().call([ocp.xVec(0),ocp.uVec(0),ocp.pVec()])
    ocp.constrain(c0,'==',0)

    # bounds
    r = 0.3
    ocp.bound('x',(-0.5,0.5))
    ocp.bound('z',(-0.5,0.5))
    ocp.bound('dx',(-5,5))
    ocp.bound('dz',(-5,5))
    ocp.bound('torque',(-50,50))
    ocp.bound('m',(0.3,0.5))
    ocp.bound('endTime',(0.01,1.5))

    # boundary conditions
    ocp.bound('x',(r,r),timestep=0)
    ocp.bound('z',(0,0),timestep=0)
    ocp.bound('x',(0,0),timestep=-1)
    ocp.bound('z',(-r*1.5,-r/2),timestep=-1)
    ocp.bound('dx',(0,0),timestep=0)
    ocp.bound('dz',(0,0),timestep=0)
    ocp.bound('dx',(0,0),timestep=-1)
    ocp.bound('dz',(0,0),timestep=-1)

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
          xOpt = numpy.array(f.input(C.NLP_X_OPT))

          self.iter = self.iter + 1
          (xup,blah) = ocp.devectorize(xOpt)
          
          po = kite_pb2.PendulumOpt()
          po.x.extend(list(xup['x']))
          po.z.extend(list(xup['z']))
          po.endTime = xup['endTime']
          po.iters = self.iter
          publisher.send_multipart(["pendulum-opt", po.SerializeToString()])
        
    # solver
    solverOptions = [ ("linear_solver","ma27")
#                    , ("derivative_test","first-order")
                    , ("expand_f",True)
                    , ("expand_g",True)
                    , ("generate_hessian",True)
#                    , ("max_iter",1000)
                    , ("tol",1e-4)
                    ]

    constraintFunOptions = [('numeric_jacobian',False)]

    # initial conditions
    ocp.guessX([r,0,0,0])
    ocp.guess('torque',0)
    ocp.guess('m',0)
    ocp.guess('endTime',0.3)

    (opt,blah) = ocp.run( ocp.lookup('endTime'),
                          solverOpts=solverOptions,
                          constraintFunOpts=constraintFunOptions,
                          callback=MyCallback() )

    # Plot the results
    plt.figure(1)
    plt.clf()
    legend = []
    for name in ocp.dae.xNames():
        legend.append(name)
        plt.plot(blah['tgrid'],opt[name])#,'--')
    for name in ocp.dae.uNames():
        legend.append(name)
        plt.plot(blah['tgrid'],opt[name]/20)#,'--')
    plt.title("pendulum swingup optimization")
    plt.xlabel('time')
    plt.legend(legend)
    plt.grid()
    plt.show()

if __name__=='__main__':
    main()
