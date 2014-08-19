#!/usr/bin/env ipython

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
import load_data
from rawe.models.arianne_conf import makeConf

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

def setupOcp(dae,conf,endTime,nk,nicp=1,deg=4):
    ocp = rawe.collocation.Coll(dae, nk=nk,nicp=nicp,deg=deg)

    print "setting up collocation..."
    ocp.setupCollocation(endTime)

    # constrain invariants
    def constrainInvariantErrs():
        R_c2b = ocp.lookup('R_c2b',timestep=0)
        rawekite.kiteutils.makeOrthonormal(ocp, R_c2b)
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0, tag=('initial c 0',None))
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0, tag=('initial cdot 0',None))
        ocp.constrain(ocp('sin_delta',timestep=0)**2 + ocp('cos_delta',timestep=0)**2,
                      '==', 1, tag=('sin**2 + cos**2 == 1',None))
    constrainInvariantErrs()

    # constrain line angle
#    for k in range(0,nk):
#        ocp.constrain(ocp.lookup('cos_line_angle',timestep=k),'>=',C.cos(45*pi/180), tag=('line angle',k))

    # constrain airspeed
    constrainCL = False
    constrainAirspeed = True
    constrainAlphaBeta = True
    def constrainAirspeed():
        for k in range(0,nk):
            for j in range(1,deg+1):
                if constrainAirspeed:
                    ocp.constrainBnds(ocp.lookup('airspeed',timestep=k,degIdx=j),
                                      (5,35), tag=('airspeed',(k,j)))
                if constrainCL:
                    ocp.constrain(ocp.lookup('cL',timestep=k,degIdx=j), '>=', 0, tag=('CL positive',nk))
                if constrainAlphaBeta:
                    ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k,degIdx=j),
                                      (-4.5,15), tag=('alpha(deg)',nk))
                    ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k,degIdx=j),
                                      (-15,15), tag=('beta(deg)',nk))

    # constrain tether force
#    constrainTetherForce(ocp)
#    realMotorConstraints(ocp)

    # bounds
    #ocp.bound('aileron', (numpy.radians(-10),numpy.radians(10)))
    #ocp.bound('elevator',(numpy.radians(-10),numpy.radians(10)))
    ocp.bound('daileron',  (numpy.radians(-200), numpy.radians(200)))
    ocp.bound('delevator', (numpy.radians(-200), numpy.radians(200)))
    #ocp.bound('ddelta',(-0.98*2*pi, 0.98*2*pi))
    ocp.bound('ddelta',(-4*pi, 4*pi))

    ocp.bound('x',(0.5,2000))
    ocp.bound('y',(-2000,2000))
    ocp.bound('z',(0, 2))
    ocp.bound('r',(1.7,1.7))
    ocp.bound('dr',(-1,1))
    ocp.bound('ddr',(-15,15))
    ocp.bound('dddr',(-15,15))

    ocp.bound('motor_torque',(-500,500))
    ocp.bound('dmotor_torque',(-100,100))

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

    return ocp


if __name__=='__main__':
    nk = 100
    tStart = 20.0 # [sec]
    tEnd   = 21.0 # [sec]
    T = tEnd - tStart
    ts = T / float(nk)

    # load the data
    (data,interval,_) = load_data.load(tStart, tEnd, ts)

    # fix some signs in the data
    data['encoder']['sin_delta'] = -data['encoder']['sin_delta']

    conf = makeConf()
    conf['useVirtualForces'] = 'random_walk'
    conf['useVirtualTorques'] = 'random_walk'
    dae = carousel_dae.makeDae(conf)

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,T,nk,deg=3)

    lineRadiusGuess = 1.7

    # trajectory for initial guess
    delta0 = numpy.arctan2(data['encoder']['sin_delta'][0], data['encoder']['cos_delta'][0])
    ddelta0 = -data['encoder']['speed_rpm'][0]*2*pi/60

    def get_steady_state_guess():
        from rawekite.carouselSteadyState import getSteadyState
        conf_ss = makeConf()
        steadyState,_ = getSteadyState(carousel_dae.makeDae(conf_ss), None, ddelta0, lineRadiusGuess, {})
        return steadyState
    ss = get_steady_state_guess()
    print ss
    ss['f1_disturbance'] = 0
    ss['f2_disturbance'] = 0
    ss['f3_disturbance'] = 0
    ss['t1_disturbance'] = 0
    ss['t2_disturbance'] = 0
    ss['t3_disturbance'] = 0

    ss['df1_disturbance'] = 0
    ss['df2_disturbance'] = 0
    ss['df3_disturbance'] = 0
    ss['dt1_disturbance'] = 0
    ss['dt2_disturbance'] = 0
    ss['dt3_disturbance'] = 0

    ocp.bound('f1_disturbance', (-100,100))
    ocp.bound('f2_disturbance', (-100,100))
    ocp.bound('f3_disturbance', (-100,100))
    ocp.bound('t1_disturbance', (-100,100))
    ocp.bound('t2_disturbance', (-100,100))
    ocp.bound('t3_disturbance', (-100,100))

    ocp.bound('df1_disturbance',(-100,100))
    ocp.bound('df2_disturbance',(-100,100))
    ocp.bound('df3_disturbance',(-100,100))
    ocp.bound('dt1_disturbance',(-100,100))
    ocp.bound('dt2_disturbance',(-100,100))
    ocp.bound('dt3_disturbance',(-100,100))

    # initial guess
    k = 0
    for nkIdx in range(ocp.nk+1):
        for nicpIdx in range(ocp.nicp):
            if nkIdx == ocp.nk and nicpIdx > 0:
                break
            for degIdx in range(ocp.deg+1):
                if nkIdx == ocp.nk and degIdx > 0:
                    break

                # path following
                delta = delta0 + ddelta0*T*(k+ocp.lagrangePoly.tau_root[degIdx])/float(ocp.nk*ocp.nicp)

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
    obj = 0

    sigmaLasElev = 0.002
    sigmaLasHoriz = 0.005
    sigmaDeltaError = 2.0*pi/180.0
    sigmaF = 100.
    sigmaT = 100.
    sigmaDF = 10.
    sigmaDT = 10.
    sigmaGyro = 3.0*pi/180.0
    sigmaAccel = 0.5
    for k,t in enumerate(interval):
        # line angle sensor
        las_h = ocp.lookup('lineAngle_hor',timestep=k)
        las_v = ocp.lookup('lineAngle_ver',timestep=k)
        obj += (data['line_angle_sensor']['angle_hor'][k] - las_h)**2/sigmaLasHoriz**2
        obj += (data['line_angle_sensor']['angle_ver'][k] - las_v)**2/sigmaLasElev**2

        # encoder
        s0 = ocp.lookup('sin_delta',timestep=k)
        c0 = ocp.lookup('cos_delta',timestep=k)
        s1 = data['encoder']['sin_delta'][k]
        c1 = data['encoder']['cos_delta'][k]
        sin_error = c0*s1 - c1*s0
        obj += sin_error**2/sigmaDeltaError**2

        # disturbances
        for name in ['f1_disturbance', 'f2_disturbance', 'f3_disturbance']:
            obj += ocp.lookup(name,timestep=k)**2/sigmaF**2
        for name in ['t1_disturbance', 't2_disturbance', 't3_disturbance']:
            obj += ocp.lookup(name,timestep=k)**2/sigmaT**2
        for name in ['df1_disturbance', 'df2_disturbance', 'df3_disturbance']:
            obj += ocp.lookup(name,timestep=k)**2/sigmaDF**2
        for name in ['dt1_disturbance', 'dt2_disturbance', 'dt3_disturbance']:
            obj += ocp.lookup(name,timestep=k)**2/sigmaDT**2

        # gyros
        for name in ['x','y','z']:
            gyroData = data['imu']['gyro_'+name][k]
            gyro = ocp.lookup('IMU_angular_velocity_'+name,timestep=k)
            obj += (gyro - gyroData)**2/sigmaGyro**2

        # accels
        for name in ['x','y','z']:
            accelData = data['imu']['accl_'+name][k]
            accel = ocp.lookup('IMU_acceleration_'+name,timestep=k)
            obj += (accel - accelData)**2/sigmaAccel**2

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

        #obj += daileron*daileron / (daileronSigma*daileronSigma)
        #obj += delevator*delevator / (delevatorSigma*delevatorSigma)
        #obj += ddr*ddr / (ddrSigma*ddrSigma)
        obj += dmotorTorque**2 / dmotorTorqueSigma**2

    # control measurements bounds
    for k in range(ocp.nk):
        # bound within 0.2 degrees
        d = 0.2*pi/180;
        ail = data['control_surfaces']['aileron'][k]
        elev = data['control_surfaces']['elevator'][k]

        ocp.bound('aileron',(ail-d,ail+d),timestep=k)
        ocp.bound('elevator',(elev-d,elev+d),timestep=k)

        # also penalize
        obj += (ocp.lookup( 'aileron',timestep=k)- ail)**2/d**2
        obj += (ocp.lookup('elevator',timestep=k)-elev)**2/d**2

    ocp.bound('aileron',(ail-d,ail+d),timestep=ocp.nk)
    ocp.bound('elevator',(-elev-d,elev+d),timestep=ocp.nk)

    obj /= float(ocp.nk)
    ocp.setObjective( obj )

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory, showAllPoints=False), 'carousel trajectory')
            ], printBoundViolation=False, printConstraintViolation=False)

    # solver
    solver = 'ipopt'
    #solver = 'snopt'
    if solver == 'ipopt':
        solverOptions = [("expand",True),
                         ("linear_solver","ma86"),
#                         ("linear_solver","ma57"),
#                         ("linear_solver","ma97"),
                         ("max_iter",1000),
                         ("tol",1e-9)]
    elif solver == 'snopt':
        solverOptions = [("expand",True)
#                        ,("_optimality_tolerance",1e-10)
#                        ,("_feasibility_tolerance",1e-10)
#                        ,("detect_linear",False)
                        ]
    else:
        raise Exception('unrecognized solver "'+solver+'"')

    print "setting up solver..."
    ocp.setupSolver( solverOpts=solverOptions,
                     callback=callback,
                     solver=solver
                   )

    xInit = None
    traj = ocp.solve(xInit=xInit)
    traj.save("data/carousel_homotopy.dat")

    # Plot the results
    def plotResults():
        interval0 = interval - interval[0]
        traj.subplot([['f1_disturbance','f2_disturbance','f3_disturbance'],
                      ['t1_disturbance','t2_disturbance','t3_disturbance']])
        traj.subplot([['df1_disturbance','df2_disturbance','df3_disturbance'],
                      ['dt1_disturbance','dt2_disturbance','dt3_disturbance']])
        #traj.subplot(['x','y','z'])
        #traj.subplot(['dx','dy','dz'])
        #traj.subplot(['r','dr','ddr'])
        #traj.subplot(['c','cdot','cddot'],title="invariants")
        #traj.subplot(['ddelta','motor_torque','dmotor_torque'],title="rotational stuff")

        # measurement errors
        traj.subplot(['sin_delta','cos_delta'])
        plt.subplot(211)
        plt.plot(interval0, data['encoder']['sin_delta'],'o')
        plt.subplot(212)
        plt.plot(interval0, data['encoder']['cos_delta'],'o')

        traj.subplot(['lineAngle_hor','lineAngle_ver'])
        plt.subplot(211)
        plt.plot(interval0, data['line_angle_sensor']['angle_hor'],'o')
        plt.subplot(212)
        plt.plot(interval0, data['line_angle_sensor']['angle_ver'],'o')

        traj.subplot(['IMU_angular_velocity_'+xyz for xyz in ['x','y','z']])
        plt.subplot(311)
        plt.plot(interval0, data['imu']['gyro_x'],'o')
        plt.subplot(312)
        plt.plot(interval0, data['imu']['gyro_y'],'o')
        plt.subplot(313)
        plt.plot(interval0, data['imu']['gyro_z'],'o')

        traj.subplot(['IMU_acceleration_'+xyz for xyz in ['x','y','z']])
        plt.subplot(311)
        plt.plot(interval0, data['imu']['accl_x'],'o')
        plt.subplot(312)
        plt.plot(interval0, data['imu']['accl_y'],'o')
        plt.subplot(313)
        plt.plot(interval0, data['imu']['accl_z'],'o')

        traj.subplot(['aileron','elevator'])
        plt.subplot(211)
        plt.plot(interval0, data['control_surfaces']['aileron'],'o')
        plt.subplot(212)
        plt.plot(interval0, data['control_surfaces']['elevator'],'o')

        #traj.subplot(['daileron','delevator'])

        traj.subplot([['airspeed'],['alpha_deg'],['beta_deg']])
        traj.subplot(['cL','cD','L_over_D'])
        #traj.subplot(['mechanical_winch_power', 'tether_tension'])
        #traj.subplot([['motor_torque'],['dmotor_torque']])
        #traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        #traj.plot(['nu'])
        plt.show()
    plotResults()
