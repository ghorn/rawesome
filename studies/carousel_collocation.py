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

import copy
import casadi as C
import matplotlib.pyplot as plt
import numpy
from numpy import pi
import pickle

import rawe
import rawekite

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
                , 3.874600000000
                , 0.0
                ])
x0=C.veccat([x0,C.sqrt(C.sumAll(x0[0:2]*x0[0:2])),0,0,0])

def setupOcp(dae,conf,nk=50,nicp=1,deg=4):
    ocp = rawe.collocation.Coll(dae, nk=nk,nicp=nicp,deg=deg)
    ocp.setupCollocation(ocp.lookup('endTime'))

    # constrain invariants
    def constrainInvariantErrs():
        dcm = ocp.lookup('dcm',timestep=0)
        rawekite.kiteutils.makeOrthonormal(ocp, dcm)
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0)
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0)
    constrainInvariantErrs()

    # constraint line angle
    for k in range(0,nk):
        ocp.constrain(ocp.lookup('cos_line_angle',timestep=k),'>=',C.cos(55*pi/180), tag=('line angle',k))

    # constrain airspeed
    def constrainAirspeedAlphaBeta():
        for k in range(0,nk):
            ocp.constrain(ocp.lookup('airspeed',timestep=k), '>=', 5)
            ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k), (-5,10))
            ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k), (-10,10))
    constrainAirspeedAlphaBeta()

    # constrain tether force
    for k in range(nk):
        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=1), '>=', 0)
        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=ocp.deg), '>=', 0)

    # make it periodic
    for name in [ "y","z",
                  "dy","dz",
                  "w1","w2","w3",
                  "ddelta",
                  "aileron","elevator",
                  "r","dr"]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1))

    # periodic attitude
#    rawekite.kiteutils.periodicEulers(ocp)
#    rawekite.kiteutils.periodicOrthonormalizedDcm(ocp)
    rawekite.kiteutils.periodicDcm(ocp)

    # bounds
    ocp.bound('aileron',(-0.04,0.04))
    ocp.bound('elevator',(-0.1,0.1))
    ocp.bound('daileron',(-2,2))
    ocp.bound('delevator',(-2,2))

    ocp.bound('x',(0.1,1000))
    ocp.bound('y',(-100,100))
    ocp.bound('z',(-0.5,7))
    ocp.bound('r',(4,4))
    ocp.bound('dr',(-10,10))
    ocp.bound('ddr',(-2.5,2.5))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-70,70))

    for w in ['w1','w2','w3']:
        ocp.bound(w,(-4*pi,4*pi))

    ocp.bound('delta',(-0.01,1.01*2*pi))
    ocp.bound('ddelta',(-pi/8,8*pi))
    ocp.bound('motor_torque',(-1000,1000))
    ocp.bound('endTime',(0.5,7.0))
    ocp.bound('w0',(10,10))

    # boundary conditions
    ocp.bound('delta',(0,0),timestep=0)
    ocp.bound('delta',(2*pi,2*pi),timestep=-1)

    # objective function
    obj = 0
    for k in range(nk):
        # control regularization
        ddr = ocp.lookup('ddr',timestep=k)
        tc = ocp.lookup('motor_torque',timestep=k)
        daileron = ocp.lookup('daileron',timestep=k)
        delevator = ocp.lookup('delevator',timestep=k)

        daileronSigma = 0.1
        delevatorSigma = 0.1
        ddrSigma = 5.0
        torqueSigma = 1.0

#        tc = tc - 390

        ailObj = daileron*daileron / (daileronSigma*daileronSigma)
        eleObj = delevator*delevator / (delevatorSigma*delevatorSigma)
        winchObj = ddr*ddr / (ddrSigma*ddrSigma)
        torqueObj = tc*tc / (torqueSigma*torqueSigma)

        obj += ailObj + eleObj + winchObj + torqueObj
    ocp.setObjective( obj/nk )

    ocp.setQuadratureDdt('quadrature_energy', 'winch_power')

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(ocp, conf, callbacks=[(rawekite.kiteTelemetry.showAllPoints,'multi-carousel')])

    # solver
    solverOptions = [ ("expand",True)
#                     ,("qp_solver",C.NLPQPSolver)
#                     ,("qp_solver_options",{'nlp_solver': C.IpoptSolver, "nlp_solver_options":{"linear_solver":"ma57"}})
                    , ("linear_solver","ma57")
                    , ("max_iter",1000)
                    , ("tol",1e-9)
#                    , ("Timeout", 1e6)
#                    , ("UserHM", True)
#                    , ("ScaleConIter",True)
#                    , ("ScaledFD",True)
#                    , ("ScaledKKT",True)
#                    , ("ScaledObj",True)
#                    , ("ScaledQP",True)
                    ]

    # initial guess
    ocp.guessX(x0)
    for k in range(0,nk+1):
        val = 2.0*pi*k/nk
        ocp.guess('delta',val,timestep=k,quiet=True)

    ocp.guess('motor_torque',0)
    ocp.guess('endTime',1.5)
    ocp.guess('daileron',0)
    ocp.guess('delevator',0)

    ocp.guess('ddr',0)
    ocp.guess('w0',10)

    ocp.setupSolver( solverOpts=solverOptions,
                     callback=callback )

    return ocp


if __name__=='__main__':
    print "reading config..."
    from highwind_carousel_conf import conf

    print "creating model..."
    dae = rawe.models.carousel(conf)
    dae.addP('endTime')

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,nk=30)

#    ocp.interpolateInitialGuess("data/carousel_opt.dat",force=True,quiet=True)

    for w0 in [10]:
        ocp.bound('w0',(w0,w0),force=True)
        traj = ocp.solve()

    print "saving optimal trajectory"
    traj.save("data/carousel_opt.dat")

    # Plot the results
    def plotResults():
        traj.plot(['x','y','z'])
        traj.plot(['aileron','elevator'],title='control surface inputs')
        traj.plot(['ddr'],title='winch accel (ddr)')
        traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed')
        traj.subplot([['alpha_deg','alphaTail_deg'],['beta_deg','betaTail_deg']])
        traj.subplot(['cL','cD','L_over_D'])
        traj.subplot(['motor_torque','motor_power'])
        traj.subplot(['winch_power','tether_tension'])
        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        plt.show()
    plotResults()
