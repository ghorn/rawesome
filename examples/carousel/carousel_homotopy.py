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
import carousel_dae
from autogen.tocarouselProto import toProto
from autogen.carousel_pb2 import Trajectory

def setupOcp(dae,conf,nk,nicp=1,deg=4):
    ocp = rawe.collocation.Coll(dae, nk=nk,nicp=nicp,deg=deg)

    print "setting up collocation..."
    ocp.setupCollocation(ocp.lookup('endTime'))

    # constrain invariants
    def constrainInvariantErrs():
        R_c2b = ocp.lookup('R_c2b',timestep=0)
        rawekite.kiteutils.makeOrthonormal(ocp, R_c2b)
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0, tag=('initial c 0',None))
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0, tag=('initial cdot 0',None))
    constrainInvariantErrs()

    # constrain line angle
    for k in range(0,nk):
        ocp.constrain(ocp.lookup('cos_line_angle',timestep=k),'>=',C.cos(55*pi/180), tag=('line angle',k))

    # constrain airspeed
    def constrainAirspeedAlphaBeta():
        for k in range(0,nk):
            for j in range(0,deg+1):
                ocp.constrainBnds(ocp.lookup('airspeed',timestep=k,degIdx=j),
                                  (10,30), tag=('airspeed',(k,j)))
                ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k,degIdx=j), (-4.5,8.5), tag=('alpha',nk))
                ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k,degIdx=j), (-6,6), tag=('beta',nk))
    constrainAirspeedAlphaBeta()

    # constrain tether force
    for k in range(nk):
        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=1), '>=', 0, tag=('tether tension',(nk,0)))
        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=ocp.deg), '>=', 0, tag=('tether tension',(nk,1)))

    # make it only go around once
    for k in range(nk):
        totalDelta = ocp.lookup('ddelta',timestep=k)*ocp.lookup('endTime')
        ocp.constrainBnds( totalDelta, (2*pi*0.9, 2*pi*1.1), tag=('total delta',k))

    # make it periodic
    for name in [ "y","z",
                  "dy","dz",
                  "w_bn_b_x","w_bn_b_y","w_bn_b_z",
                  "r","dr",'ddr',
                  'ddelta',
                  'aileron','elevator','rudder','flaps',
                  'cos_delta','sin_delta','motor_torque'
                  ]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1), tag=('periodic '+name,None))

    # periodic attitude
    rawekite.kiteutils.periodicDcm(ocp)

    # bounds
    ocp.bound('aileron', (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('elevator',(numpy.radians(-10),numpy.radians(10)))
    ocp.bound('rudder',  (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('flaps',  (numpy.radians(0),numpy.radians(0)))
    # can't bound flaps==0 AND have periodic flaps at the same time
    # bounding flaps (-1,1) at timestep 0 doesn't really free them, but satisfies LICQ
    ocp.bound('flaps', (-1,1),timestep=0,quiet=True)
    ocp.bound('daileron',  (numpy.radians(-20), numpy.radians(20)))
    ocp.bound('delevator', (numpy.radians(-20), numpy.radians(20)))
    ocp.bound('drudder',   (numpy.radians(-20), numpy.radians(20)))
    ocp.bound('dflaps',    (numpy.radians(-20), numpy.radians(20)))
    ocp.bound('ddelta',(0.0,4*pi))

    ocp.bound('x',(-2000,2000))
    ocp.bound('y',(-2000,2000))
    ocp.bound('z',(-2,1.5))
    ocp.bound('r',(4,10))
    ocp.bound('dr',(-1,1))
    ocp.bound('ddr',(-15,15))
    ocp.bound('dddr',(-15,15))

    ocp.bound('motor_torque',(-500,500))
    ocp.bound('dmotor_torque',(-1,1))

    ocp.bound('cos_delta',(-1.5,1.5))
    ocp.bound('sin_delta',(-1.5,1.5))

    ocp.bound('sin_delta', (0,0), timestep=0, quiet=True)
    ocp.bound('cos_delta', (0.5,1.5), timestep=0, quiet=True)

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-70,70))

    for w in ['w_bn_b_x',
              'w_bn_b_y',
              'w_bn_b_z']:
        ocp.bound(w,(-4*pi,4*pi))

    ocp.bound('endTime',(1.0,2.0))
    ocp.guess('endTime',1.5)
    ocp.bound('w0',(0,0))
#    ocp.bound('w0',(10,10))

    # boundary conditions
    ocp.bound('y',(0,0),timestep=0,quiet=True)

    return ocp


if __name__=='__main__':
    from rawe.models.betty_conf import makeConf
    conf = makeConf()
    conf['runHomotopy'] = True
    conf['minAltitude'] = 0.5
    nk = 40
    dae = carousel_dae.makeDae(conf)
    dae.addP('endTime')

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,nk)

    lineRadiusGuess = 4.0
    nTurns = 1

    # trajectory for homotopy
    homotopyTraj = {'sin_delta':[],'cos_delta':[]}
    k = 0
    for nkIdx in range(ocp.nk+1):
        for nicpIdx in range(ocp.nicp):
            if nkIdx == ocp.nk and nicpIdx > 0:
                break
            for degIdx in range(ocp.deg+1):
                if nkIdx == ocp.nk and degIdx > 0:
                    break

                # path following
                delta = 2*pi*(k+ocp.lagrangePoly.tau_root[degIdx])/float(ocp.nk*ocp.nicp)

                if nicpIdx == 0 and degIdx == 0:
                    homotopyTraj['sin_delta'].append(numpy.sin(delta))
                    homotopyTraj['cos_delta'].append(numpy.cos(delta))

                ocp.guess('sin_delta',numpy.sin(delta),timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('cos_delta',numpy.cos(delta),timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
            k += 1

    ocp.guess('w_bn_b_z', -2.0*pi/ocp._guess.lookup('endTime'))
    ocp.guess('ddelta', nTurns*2*pi/(ocp._guess.lookup('endTime')))
    ocp.guess('e11',0); ocp.guess('e12',1); ocp.guess('e13',0)
    ocp.guess('e21',0); ocp.guess('e22',0); ocp.guess('e23',-1)
    ocp.guess('e31',-1); ocp.guess('e32',0); ocp.guess('e33',0)

    # objective function
    obj = -1e6*ocp.lookup('gamma_homotopy')
#    mean_x = numpy.mean(homotopyTraj['x'])
#    mean_y = numpy.mean(homotopyTraj['y'])
#    mean_z = numpy.mean(homotopyTraj['z'])
#    for k in range(ocp.nk+1):
#        x = ocp.lookup('x',timestep=k)
#        y = ocp.lookup('y',timestep=k)
#        z = ocp.lookup('z',timestep=k)
#        obj += ((x-mean_x)**2 + (y-mean_y)**2 + (z-mean_z)**2 - circleRadiusGuess**2)**2

    for k in range(ocp.nk+1):
        obj += (homotopyTraj['sin_delta'][k] - ocp.lookup('sin_delta',timestep=k))**2
        obj += (homotopyTraj['cos_delta'][k] - ocp.lookup('cos_delta',timestep=k))**2
    ocp.setQuadratureDdt('mechanical_energy', 'mechanical_winch_power')
    ocp.setQuadratureDdt('electrical_energy', 'electrical_winch_power')

    # control regularization
    for k in range(ocp.nk):
        ddr = ocp.lookup('ddr',timestep=k)
        daileron = ocp.lookup('daileron',timestep=k)
        delevator = ocp.lookup('delevator',timestep=k)
        drudder = ocp.lookup('drudder',timestep=k)
        dflaps = ocp.lookup('dflaps',timestep=k)
        dmotorTorque = ocp.lookup('dmotor_torque',timestep=k)

        daileronSigma = 0.01
        delevatorSigma = 0.1
        drudderSigma = 0.1
        dflapsSigma = 0.1
        ddrSigma = 0.1
        dmotorTorqueSigma = 0.1

        ailObj = daileron*daileron / (daileronSigma*daileronSigma)
        eleObj = delevator*delevator / (delevatorSigma*delevatorSigma)
        rudObj = drudder*drudder / (drudderSigma*drudderSigma)
        flapsObj = dflaps*dflaps / (dflapsSigma*dflapsSigma)
        winchObj = ddr*ddr / (ddrSigma*ddrSigma)
        motorTorqueObj = dmotorTorque*dmotorTorque / (dmotorTorqueSigma*dmotorTorqueSigma)

        obj += 1e-2*(ailObj + eleObj + rudObj + flapsObj + winchObj + motorTorqueObj)/float(ocp.nk)

    # homotopy forces/torques regularization
    homoReg = 0
    for k in range(ocp.nk):
        for nicpIdx in range(ocp.nicp):
            for degIdx in range(1,ocp.deg+1):
                homoReg += ocp.lookup('f1_homotopy',timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)**2
                homoReg += ocp.lookup('f2_homotopy',timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)**2
                homoReg += ocp.lookup('f3_homotopy',timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)**2
                homoReg += ocp.lookup('t1_homotopy',timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)**2
                homoReg += ocp.lookup('t2_homotopy',timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)**2
                homoReg += ocp.lookup('t3_homotopy',timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)**2
    obj += 1e-2*homoReg/float(ocp.nk*ocp.nicp*ocp.deg)

    ocp.setObjective( obj )

    # initial guesses
    ocp.guess('w0',0)
    ocp.guess('r',lineRadiusGuess)
    ocp.guess('x',lineRadiusGuess)

    for name in ['w_bn_b_x','w_bn_b_y','dr','ddr','dddr','aileron','rudder','flaps','motor_torque',
                 'elevator','daileron','delevator','drudder','dflaps','dmotor_torque',
                 'y','z','dx','dy','dz']:
        ocp.guess(name,0)

    ocp.guess('gamma_homotopy',0)

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory, showAllPoints=False), 'carousel trajectory')
            ], printBoundViolation=True, printConstraintViolation=True)

    # solver
    solverOptions = [("linear_solver","ma27"),
                     ("max_iter",1000),
                     ("expand",True),
                     ("tol",1e-10)]

    print "setting up solver..."
    ocp.setupSolver( solverOpts=solverOptions,
                     callback=callback )

    xInit = None
    ocp.bound('gamma_homotopy',(1e-4,1e-4),force=True)
    traj = ocp.solve(xInit=xInit)

    ocp.bound('gamma_homotopy',(0,1),force=True)
    traj = ocp.solve(xInit=traj.getDvs())

    ocp.bound('gamma_homotopy',(1,1),force=True)
#    ocp.bound('endTime',(3.5,6.0),force=True)
    traj = ocp.solve(xInit=traj.getDvs())

    traj.save("data/crosswind_homotopy.dat")

    def printBoundsFeedback():
        # bounds feedback
        xOpt = traj.dvMap.vectorize()
        lbx = ocp.solver.input('lbx')
        ubx = ocp.solver.input('ubx')
        ocp._bounds.printBoundsFeedback(xOpt,lbx,ubx,reportThreshold=0)
    printBoundsFeedback()

    # Plot the results
    def plotResults():
        traj.subplot(['f1_homotopy','f2_homotopy','f3_homotopy'])
        traj.subplot(['t1_homotopy','t2_homotopy','t3_homotopy'])
        traj.subplot(['x','y','z'])
        traj.subplot(['dx','dy','dz'])
        traj.subplot([['aileron','elevator'],['daileron','delevator']],title='control surfaces')
        traj.subplot([['rudder','flaps'],['drudder','dflaps']],title='control surfaces')
        traj.subplot(['r','dr','ddr'])
        traj.subplot(['wind_at_altitude','dr','dx'])
        traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed')
        traj.subplot([['alpha_deg'],['beta_deg']])
        traj.subplot(['cL','cD','L_over_D'])
        traj.subplot(['mechanical_winch_power', 'tether_tension'])
        traj.subplot([['motor_torque'],['dmotor_torque']])
        traj.subplot(['w_bn_b_x','w_bn_b_y','w_bn_b_z'])
        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        traj.plot(['nu'])
        plt.show()
    plotResults()
