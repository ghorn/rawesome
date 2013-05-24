import matplotlib.pyplot as plt

from rawe.dae import Dae
from rawe.collocation import Coll

if __name__ == "__main__":
    ######## make the Dae #######
    dae = Dae()

    [pos,vel,mass] = dae.addX( ["pos","vel","mass"] )
    thrust = dae.addU( "thrust" )
    
    # some extra outputs for the dae model
    dae['posvel'] = pos*vel
    dae['velvel'] = vel*vel

    # specify the ode residual
    dae.setResidual([dae.ddt('pos') - vel,
                     dae.ddt('vel') - thrust/mass,
                     dae.ddt('mass') - 0.1*thrust*thrust])

    ######## make the collocation scheme ########
    N = 100
    ocp = Coll(dae, nk=N, nicp=1, collPoly="RADAU", deg=4)

    endTime = 5
    ocp.setupCollocation( endTime )

    # bounds
    ocp.bound('thrust',(-1.3,0.9))
    ocp.bound('pos',(-10,10))
    ocp.bound('vel',(-10,10))
    ocp.bound('mass',(0.001,1000))

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
    posvel4_2 = ocp('posvel',timestep=4, degIdx=2)

    # add nonlinear constraint
    ocp.constrain(thrust4, '<=', posvel4_2**2)

    # fix objective and setup solver
    obj = sum([ocp('thrust',timestep=k)**2
               for k in range(ocp.nk)])
    ocp.setObjective(obj)
#    ocp.setObjective(ocp.lookup('integral vel*vel',timestep=-1))

    solverOptions = [ ("tol",1e-9) ]
    ocp.setupSolver(solverOpts=solverOptions)

#    ocp.interpolateInitialGuess("data/rocket_opt.dat",force=True,quiet=True,numLoops=1)
    traj = ocp.solve()

#    maxErr = 0.0
#    for k in range(ocp.nk):
#        for d in range(ocp.deg+1):
#            n1 = traj.lookup('integral vel*vel',timestep=k,nicpIdx=0,degIdx=d)
#            n2 = traj.lookup('integral vel*vel2',timestep=k,nicpIdx=0,degIdx=d)
#            err = 100.0*(n1-n2)/(n1+1e-12)
#            if abs(err) > abs(maxErr):
#                maxErr = float(err)
#    n1 = traj.lookup('integral vel*vel',timestep=-1)
#    n2 = traj.lookup('integral vel*vel2',timestep=-1)
#    err = 100.0*(n1-n2)/(n1+1e-12)
#    if abs(err) > abs(maxErr):
#        maxErr = float(err)
#    print "maxErr: "+str(err)+" %"


    print "final position: "+str(traj.lookup('pos',-1))
    
    # save trajectory
    traj.save("data/rocket_opt.dat")

    # plot results
#    traj.plot('pos')
    traj.plot(['pos','vel'])
    traj.plot('thrust')
#    traj.subplot([['pos','vel'],['thrust']])
#    traj.plot('pos*vel')
#    traj.subplot(['integral vel*vel','integral vel*vel2'])
#    traj.plot(['integral vel*vel','integral vel*vel2'])
    plt.show()
