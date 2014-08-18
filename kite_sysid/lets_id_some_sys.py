#!/usr/bin/env python

# Copyright 2012-2014 Greg Horn
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
import scipy.io
from numpy import pi
import pickle

import rawe
import rawekite
from autogen.tocarousel_sysidProto import toProto
from autogen.carousel_sysid_pb2 import Trajectory
import carousel_dae

#def constrainTetherForce(ocp):
#    for k in range(ocp.nk):
#        for j in range(1,ocp.deg+1): #[1]:
#            ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=j), '>=', 0, tag=('tether tension positive',k))
#
#def realMotorConstraints(ocp):
#    safety_factor = 0.5
#    for k in range(nk):
##        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=1),       '<=', 150, tag=('motor torque',k))
##        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=ocp.deg), '<=', 150, tag=('motor torque',k))
#        ocp.constrainBnds( ocp.lookup('torque',timestep=k,degIdx=1),
#                           (-78*safety_factor,78*safety_factor), tag=('motor torque',k))
#        ocp.constrainBnds( ocp.lookup('torque',timestep=k,degIdx=ocp.deg),
#                           (-78*safety_factor,78*safety_factor), tag=('motor torque',k))
#
#        ocp.constrain( ocp.lookup('rpm',timestep=k),       '<=', 1500*safety_factor, tag=('rpm',k))
#        ocp.constrain( -1500*safety_factor, '<=', ocp.lookup('rpm',timestep=k),       tag=('rpm',k))

def setupOcp(dae,conf,nk,nicp=1,deg=4):
    ocp = rawe.collocation.Coll(dae, nk=nk,nicp=nicp,deg=deg)

    print "setting up collocation..."
    ocp.setupCollocation(ocp.lookup('endTime'))

#    for k in range(ocp.nk):
#        for j in range(1,ocp.deg+1): #[1]:
#            ocp.constrain( ocp('x',timestep=k,degIdx=j), '>=', ocp('y',timestep=k,degIdx=j), tag=('x>y',k))
#            ocp.constrain( ocp('x',timestep=k,degIdx=j), '>=', -ocp('y',timestep=k,degIdx=j), tag=('x>-y',k))

    # constrain invariants
    def constrainInvariantErrs():
        R_c2b = ocp.lookup('R_c2b',timestep=0)
        rawekite.kiteutils.makeOrthonormal(ocp, R_c2b)
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0, tag=('initial c 0',None))
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0, tag=('initial cdot 0',None))
        ocp.constrain(ocp('sin_delta',timestep=0)**2 + ocp('cos_delta',timestep=0)**2,
                      '==', 1, tag=('sin**2 + cos**2 == 1',None))
    constrainInvariantErrs()

#    # constrain line angle
#    for k in range(0,nk):
#        ocp.constrain(ocp.lookup('cos_line_angle',timestep=k),'>=',C.cos(45*pi/180), tag=('line angle',k))

    # constrain airspeed
    def constrainAirspeedAlphaBeta():
        for k in range(0,nk):
            for j in range(1,deg+1):
                ocp.constrainBnds(ocp.lookup('airspeed',timestep=k,degIdx=j),
                                  (5,35), tag=('airspeed',(k,j)))
                #ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k,degIdx=j), (-4.5,7.5), tag=('alpha(deg)',nk))
                #ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k,degIdx=j), (-6,6), tag=('beta(deg)',nk))
                ocp.constrain(ocp.lookup('cL',timestep=k,degIdx=j), '>=', 0, tag=('CL positive',nk))
    constrainAirspeedAlphaBeta()

    # constrain tether force
#    constrainTetherForce(ocp)
#    realMotorConstraints(ocp)

    # make it periodic
    for name in [ "y","z",
                  "dy","dz",
                  "w_bn_b_x","w_bn_b_y","w_bn_b_z",
                  "r","dr",'ddr',
                  'ddelta',
                  'aileron','elevator',
                  'sin_delta','motor_torque'
                  ]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1), tag=('periodic '+name,None))

    # periodic attitude
    rawekite.kiteutils.periodicDcm(ocp)

    # bounds
    ocp.bound('aileron', (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('elevator',(numpy.radians(-10),numpy.radians(10)))
    ocp.bound('daileron',  (numpy.radians(-20), numpy.radians(20)))
    ocp.bound('delevator', (numpy.radians(-20), numpy.radians(20)))
    #ocp.bound('ddelta',(-0.98*2*pi, 0.98*2*pi))
    ocp.bound('ddelta',(-4*pi, 4*pi))

    ocp.bound('x',(0.5,2000))
    ocp.bound('y',(-2000,2000))
    ocp.bound('z',(0, 2))
    ocp.bound('r',(1.6,1.8))
    # can't bound r AND have periodic r at the same time
    ocp.bound('r',(0,100),timestep=-1)
    ocp.bound('dr',(-1,1))
    ocp.bound('ddr',(-15,15))
    ocp.bound('dddr',(-15,15))

    ocp.bound('motor_torque',(-500,500))
    ocp.bound('dmotor_torque',(-1,1))

    ocp.bound('cos_delta',(-1.5,1.5))
    ocp.bound('sin_delta',(-1.5,1.5))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-70,70))

    for w in ['w_bn_b_x',
              'w_bn_b_y',
              'w_bn_b_z']:
        ocp.bound(w,(-4*pi,4*pi))

    ocp.bound('endTime',(1.2, 1.2))
    ocp.guess('endTime',1.2)

    # boundary conditions
    ocp.bound('sin_delta', (0,0), timestep=0, quiet=True)
    ocp.bound('cos_delta', (0.5,1.5), timestep=0, quiet=True)
    ocp.bound('cos_delta', (0.5,1.5), timestep=-1, quiet=True)

    return ocp


if __name__=='__main__':
    datFolder = 'data/smc_20140430_232938_dmhe_testing'
    data = {}
    for name in ['cable_length','control_surfaces','encoder',
                 'imu','led_tracker','line_angle_sensor']:
        mat = scipy.io.loadmat(datFolder + '/' + name + '.mat')

    from rawe.models.arianne_conf import makeConf
    conf = makeConf()
    conf['runHomotopy'] = True
    nk = 35
    #nk = 70
    dae = carousel_dae.makeDae(conf)
    dae.addP('endTime')

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,nk)

    lineRadiusGuess = 1.7
    nTurns = 1

    # trajectory for homotopy
    def get_steady_state_guess():
        ddelta = -nTurns*2*pi/(ocp._guess.lookup('endTime'))
        from rawekite.carouselSteadyState import getSteadyState
        conf_ss = makeConf()
        steadyState,_ = getSteadyState(carousel_dae.makeDae(conf_ss), None, ddelta, lineRadiusGuess, {})
        return steadyState
    ss = get_steady_state_guess()
    print ss

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
                delta = -2*pi*(k+ocp.lagrangePoly.tau_root[degIdx])/float(ocp.nk*ocp.nicp)

                if nicpIdx == 0 and degIdx == 0:
                    homotopyTraj['sin_delta'].append(numpy.sin(delta))
                    homotopyTraj['cos_delta'].append(numpy.cos(delta))

                ocp.guess('sin_delta',numpy.sin(delta),timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('cos_delta',numpy.cos(delta),timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                for name in dae.xNames():
                    if name not in ['sin_delta','cos_delta']:
                        ocp.guess(name,ss[name],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                if degIdx > 0:
                    for name in dae.zNames():
                        if name not in ['f1_homotopy', 'f2_homotopy', 'f3_homotopy',
                                        't1_homotopy', 't2_homotopy', 't3_homotopy']:
                            ocp.guess(name,ss[name],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
            k += 1
    for nkIdx in range(ocp.nk):
        for name in dae.uNames():
            ocp.guess(name,ss[name],timestep=nkIdx)

    # objective function
    obj = -1e6*ocp.lookup('gamma_homotopy')

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
        dmotorTorque = ocp.lookup('dmotor_torque',timestep=k)

        daileronSigma = 0.01
        delevatorSigma = 0.1
        ddrSigma = 0.1
        dmotorTorqueSigma = 1.0

        ailObj = daileron*daileron / (daileronSigma*daileronSigma)
        eleObj = delevator*delevator / (delevatorSigma*delevatorSigma)
        winchObj = ddr*ddr / (ddrSigma*ddrSigma)
        motorTorqueObj = dmotorTorque*dmotorTorque / (dmotorTorqueSigma*dmotorTorqueSigma)

        obj += 1e-2*(ailObj + eleObj + winchObj + motorTorqueObj)/float(ocp.nk)

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
    ocp.guess('gamma_homotopy',0)

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory, showAllPoints=False), 'carousel trajectory')
            ], printBoundViolation=False, printConstraintViolation=False)

    # solver
    solverOptions = [("expand",True),
                     ("linear_solver","ma86"),
                     ("max_iter",1000),
                     ("tol",1e-10)
                     #("_optimality_tolerance",1e-10),
                     #("_feasibility_tolerance",1e-13),
                     #("_detect_linear",False),
                     ]

    print "setting up solver..."
    ocp.setupSolver( solverOpts=solverOptions,
                     callback=callback,
                     solver='ipopt'
                     #solver='snopt'
                   )

    xInit = None
#    ocp.bound('gamma_homotopy',(1e-4,1e-4),force=True)
    ocp.bound('gamma_homotopy',(1,1),force=True)
    traj = ocp.solve(xInit=xInit)

#    ocp.bound('gamma_homotopy',(0,1),force=True)
#    traj = ocp.solve(xInit=traj.getDvs())

#    ocp.bound('gamma_homotopy',(1,1),force=True)
##    ocp.bound('endTime',(3.5,6.0),force=True)
#    traj = ocp.solve(xInit=traj.getDvs())

    traj.save("data/carousel_homotopy.dat")

    # Plot the results
    def plotResults():
        #traj.subplot(['f1_homotopy','f2_homotopy','f3_homotopy'])
        #traj.subplot(['t1_homotopy','t2_homotopy','t3_homotopy'])
        #traj.subplot(['x','y','z'])
        #traj.subplot(['dx','dy','dz'])
        #traj.subplot([['aileron','elevator'],['daileron','delevator']],title='control surfaces')
        #traj.subplot(['r','dr','ddr'])
        #traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed')
        traj.subplot([['alpha_deg'],['beta_deg']])
        traj.subplot(['cL','cD','L_over_D'])
        #traj.subplot(['mechanical_winch_power', 'tether_tension'])
        #traj.subplot([['motor_torque'],['dmotor_torque']])
        #traj.subplot(['w_bn_b_x','w_bn_b_y','w_bn_b_z'])
        #traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        #traj.plot(['nu'])
        plt.show()
#    plotResults()
