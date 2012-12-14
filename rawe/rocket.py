from models.dae import Dae
from collocation import Coll
from trajectory import Trajectory
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # make the Dae
    dae = Dae()

    [pos,vel,mass] = dae.addX( ["pos","vel","mass"] )
    [zdummy] = dae.addZ( ["zdummy"] )
    thrust = dae.addU( "thrust" )
    
    dae['pos*vel'] = pos*vel

    dae.setOdeRes([dae.ddt('pos') - vel,
                   dae.ddt('vel') - (thrust - 0.05*vel*vel)/mass,
                   dae.ddt('mass') - -0.1*thrust*thrust])

    dae.setAlgRes([zdummy])

    # make the collocation scheme
    N = 100
    ocp = Coll(dae, nk=N, nicp=1, collPoly="RADAU", deg=4)

    endTime = 5
    ocp.setupCollocation( endTime )

    # bounds
    ocp.bound('thrust',(-0.9,0.8))
    ocp.bound('pos',(-10,10))
    ocp.bound('vel',(-10,10))
    ocp.bound('mass',(0,1000))

    # boundary conditions
    ocp.bound('pos',(0,0),timestep=0)
    ocp.bound('pos',(5,5),timestep=-1)

    ocp.bound('vel',(0,0),timestep=0)
    ocp.bound('vel',(0,0),timestep=-1)

    ocp.bound('mass',(1,1),timestep=0)

    ocp.guess("pos",0)
    ocp.guess("vel",0)
    ocp.guess("mass",1)
    ocp.guess("thrust",0)
    

    # lookup states/actions/outputs/params
    thrust4 = ocp.lookup('thrust',timestep=4)
    thrust4 = ocp('thrust',timestep=4)
    
    # can specify index of collocation point
#    posvel4_2 = ocp('pos*vel',timestep=4, degIdx=2)

    # add nonlinear constraint
#    ocp.constrain(thrust4, '<=', posvel4_2**2)

    # fix objective and setup solver
    obj = sum([ocp('thrust',timestep=k)**2
               for k in range(ocp.nk)])
    ocp.setObjective(obj)
    
    solverOptions = [ ("expand_f",True)
                    , ("expand_g",True)
                    , ("generate_hessian",True)
                    , ("tol",1e-9)
                    ]
    ocp.setupSolver(solverOpts=solverOptions)

    opt = ocp.solve()

    # make trajectory
    traj = Trajectory(ocp,dvs=opt.vec)

    print "final position: "+str(opt.lookup('pos',-1))
    
    # save trajectory
    traj.save("data/rocket_opt.dat")

    # plot results
    traj.plot('pos')
    traj.plot(['pos','vel'])
    traj.subplot([['pos','vel'],['thrust']])
    plt.show()
