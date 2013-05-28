import matplotlib.pyplot as plt
import zmq
import numpy
from numpy import pi
import copy

import casadi as C

import rawe

def main():
    nk = 50

    print "creating model"
    dae = rawe.models.pendulum2()
    dae.addP('endTime')

    print "setting up OCP"
    ocp = rawe.collocation.Coll(dae, nk=nk,nicp=1,deg=4, collPoly='RADAU')
    print "setting up collocation"
    ocp.setupCollocation( ocp.lookup('endTime') )
    
    # constrain invariants
    ocp.constrain(ocp.lookup('c',timestep=0),'==',0)
    ocp.constrain(ocp.lookup('cdot',timestep=0),'==',0)

    # bounds
    r = 0.3
    ocp.bound('x',(-2*r,2*r))
    ocp.bound('z',(-2*r,0.01*r))
    ocp.bound('dx',(-5,5))
    ocp.bound('dz',(-5,5))
    ocp.bound('torque',(-50,50))
    ocp.bound('m',(0.3,0.3))
    ocp.bound('endTime',(0.1,3.5))

    # boundary conditions
    ocp.bound('x',(r,r),timestep=0)
    ocp.bound('z',(0,0),timestep=0)
    ocp.bound('x',(0,0),timestep=-1)
    ocp.bound('z',(-10*r,0.01*r),timestep=-1)
    ocp.bound('dx',(0,0),timestep=0)
    ocp.bound('dz',(0,0),timestep=0)
    ocp.bound('dx',(0,0),timestep=-1)
    ocp.bound('dz',(-0.5,0.5),timestep=-1)

    # make the solver
    obj = 0
    for k in range(ocp.nk):
        t = ocp.lookup('torque',timestep=k)
        obj += t*t

    ocp.setObjective(ocp.lookup('endTime') + 1e-6*obj/float(nk))
#    ocp.setObjective(1e-6*obj/float(nk))
    
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
            traj = rawe.collocation.trajectory.Trajectory(ocp,xOpt)
            
            po = rawe.kite_pb2.PendulumOpt()
            po.x.extend(list([traj.lookup('x',timestep=k) for k in range(ocp.nk+1)]))
            po.z.extend(list([traj.lookup('z',timestep=k) for k in range(ocp.nk+1)]))
            po.messages.append('endTime: %.3f'% traj.lookup('endTime'))
            po.messages.append('mass: %.3f'% traj.lookup('m'))
            po.messages.append('iters: %d'  % self.iter)
            publisher.send_multipart(["pendulum-opt", po.SerializeToString()])
        
    # solver
    solverOptions = [ ("linear_solver","ma57")
#                    , ("derivative_test","first-order")
                    , ("expand_f",True)
                    , ("expand_g",True)
                    , ("generate_hessian",True)
                    , ("max_iter",10000)
                    , ("tol",1e-8)
                    ]
#    solverOptions = [ ("Timeout", 1e6),
#                      ("UserHM", True)]
#                      ("ScaleConIter",True),
#                      ("ScaledFD",True),
#                      ("ScaledKKT",True),
#                      ("ScaledObj",True),
#                      ("ScaledQP",True)]
    
    constraintFunOptions = [('numeric_jacobian',False)]

    # initial conditions
    endTime = 0.3
    xOld = r
    zOld = 0
    dt0 = endTime/nk
    for k in range(nk+1):
        theta = float(k)/nk*C.pi/8.0
        x =  r*C.cos(theta)
        z = -r*C.sin(theta)
        
        ocp.guess('x',x,timestep=k)
        ocp.guess('z',z,timestep=k)
        ocp.guess('dx',(x-xOld)/dt0,timestep=k)
        ocp.guess('dz',(z-zOld)/dt0,timestep=k)
        xOld = x
        zOld = z
    for k in range(ocp.nk):
        if k < ocp.nk:
            ocp.guess('torque',0,timestep=k)
        else:
            ocp.guess('torque',-0,timestep=k)
    ocp.guess('m',0.3)
    ocp.guess('endTime',endTime)

##    ocp.guess('dx',0)
##    ocp.guess('dz',0)
#    ocp.guessX([r,0,0,0])
#    ocp.guessX([0,-r,0,0])
#    ocp.guess('torque',0)

#    ocp.interpolateInitialGuess("data/pendulum_opt.dat",force=True,quiet=True)

    print "setting up solver"
    ocp.setupSolver( solverOpts=solverOptions,
                     constraintFunOpts=constraintFunOptions,
                     callback=MyCallback() )
    
    print "solving"
    traj = ocp.solve()
    print "endTime: "+str(traj.lookup('endTime'))
    print "mass: "+str(traj.lookup('m'))

    print "saving optimal trajectory"
    traj.save("data/pendulum_opt.dat")

    # Plot the results
    traj.subplot([['x','z'],['dx','dz']])
    traj.plot('torque')
    traj.plot('tau')
    plt.show()

if __name__=='__main__':
    main()
