import matplotlib.pyplot as plt

from rawe.dae import Dae
from rawe.collocation import Coll

if __name__ == "__main__":
    ######## make the Dae #######
    dae = Dae()

    [pos,vel] = dae.addX( ["pos","vel"] )
    dae.addX('abspos')
    force = dae.addU( "force" )
    
    # some extra outputs for the dae model
    dae['force/pos'] = force/pos

    # specify the ode residual
    mass = 1.0
    dae.setResidual([dae.ddt('pos') - vel,
                     dae.ddt('vel') - (force/mass - 3.0*pos - 0.2*vel)]
                    )

    ######## make the collocation scheme ########
    N = 100
    ocp = Coll(dae, nk=N, nicp=1, collPoly="RADAU", deg=4)

    endTime = 5
    ocp.setupCollocation( endTime )

    # bounds
    ocp.bound('abspos',(0,10))
    ocp.bound('pos',(-1000,1000))
    ocp.bound('vel',(-100,100))
    ocp.bound('force',(-30,30))

    # boundary conditions
    ocp.bound('pos',(5,5),timestep=0)
    ocp.bound('pos',(0,0),timestep=-1)

    ocp.bound('vel',(0,0),timestep=0)
    ocp.bound('vel',(0,0),timestep=-1)

    ocp.guess("abspos",0)
    ocp.guess("pos",0)
    ocp.guess("vel",0)
    ocp.guess("force",0)
    
    # add slacks
    for k in range(N):
        for j in range(ocp.deg+1):
            pos =    ocp('pos',   timestep=k, degIdx=j)
            abspos = ocp('abspos',timestep=k, degIdx=j)
            
            ocp.constrain(    pos, '<=', abspos)
            ocp.constrain(-abspos, '<=',    pos)
    pos =    ocp('pos',   timestep=N, degIdx=0)
    abspos = ocp('abspos',timestep=N, degIdx=0)
    ocp.constrain(    pos, '<=', abspos)
    ocp.constrain(-abspos, '<=',    pos)

    # fix objective and setup solver
    # L1 parts
    obj = 0
    for k in range(N):
        for j in range(ocp.deg+1):
            obj += ocp('abspos',timestep=k,degIdx=j)
    obj += ocp('abspos',timestep=N,degIdx=0)

    # L2 parts
    for k in range(N):
        obj += ocp('force',timestep=k)**2

    ocp.setObjective(obj)

    solverOptions = [ ("expand_f",True)
                    , ("expand_g",True)
                    , ("generate_hessian",True)
                    , ("tol",1e-9)
                    ]
    ocp.setupSolver(solverOpts=solverOptions)

    traj = ocp.solve()

    # plot results
    traj.subplot([['pos','vel'],['force'],['force/pos']])
    traj.plot(['pos','abspos'])
    plt.show()
