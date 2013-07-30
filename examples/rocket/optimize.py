# Copyright 2012-2013 Greg Horn
#
# This file is part of rawesome.
#
# rawesome is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rawesome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib.pyplot as plt

import rawe
from rawe.collocation import Coll
import rocket_dae
from autogen.torocketProto import toProto
from autogen.rocket_pb2 import Trajectory

if __name__ == "__main__":
    dae = rocket_dae.makeDae()

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

    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory), 'rocket trajectory')
        ])

    solverOptions = [ ("tol",1e-9) ]
    ocp.setupSolver(solverOpts=solverOptions,
                    callback=callback)

#    ocp.interpolateInitialGuess("data/rocket_opt.dat",force=True,quiet=True,numLoops=1)
    traj = ocp.solve()


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
