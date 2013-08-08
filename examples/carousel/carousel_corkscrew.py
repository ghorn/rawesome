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
                                  (10,45), tag=('airspeed',(k,j)))
                ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k,degIdx=j), (-4.5,8.5), tag=('alpha',(k,j)))
                ocp.constrain(ocp.lookup('cL',timestep=k,degIdx=j), '>=', -0.1, tag=('CL > -0.1',(k,j)))
                ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k,degIdx=j), (-6,6), tag=('beta',(k,j)))
                x = ocp('x', timestep=k,degIdx=j)
                y = ocp('y', timestep=k,degIdx=j)
                z = ocp('z', timestep=k,degIdx=j)
                ocp.constrain(C.sqrt(x**2 + y**2), '>=', -z, tag=('azimuth not too high',(k,j)))
    constrainAirspeedAlphaBeta()

    # constrain tether force
    for k in range(nk):
        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=1), '>=', 0, tag=('tether tension',(nk,0)))
        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=ocp.deg), '>=', 0, tag=('tether tension',(nk,1)))

    # bounds
    ocp.bound('aileron', (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('elevator',(numpy.radians(-10),numpy.radians(10)))
    ocp.bound('rudder',  (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('flaps',  (numpy.radians(0),numpy.radians(0)))
    ocp.bound('daileron',  (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('delevator', (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('drudder',   (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('dflaps',    (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('ddelta',(-4*pi,4*pi))

    ocp.bound('x',(-2000,2000))
    ocp.bound('y',(-2000,2000))
    ocp.bound('z',(-50,1.5))
    ocp.bound('r',(4,300))
    ocp.bound('dr',(-10,10))
    ocp.bound('ddr',(-15,15))
    ocp.bound('dddr',(-50,50))

    ocp.bound('motor_torque',(-500,500))
    ocp.bound('dmotor_torque',(-100,100))

    ocp.bound('cos_delta',(-1.5,1.5))
    ocp.bound('sin_delta',(-1.5,1.5))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-1000,1000))

    for w in ['w_bn_b_x',
              'w_bn_b_y',
              'w_bn_b_z']:
        ocp.bound(w,(-4*pi,4*pi))

    ocp.bound('endTime',(1,50))
#    ocp.guess('endTime',3.5)
    ocp.bound('w0',(4,4))
#    ocp.bound('w0',(10,10))

    # boundary conditions
    ocp.bound('cos_delta', (0.5,1.5), timestep=0, quiet=True)

    # boundary conditions
    ocp.bound('sin_delta', (0,0), timestep=-1, quiet=True)
    ocp.bound('cos_delta', (0.5,1.5), timestep=-1, quiet=True)
    ocp.bound('ddelta', (0,0), timestep=-1, quiet=True)

    f = open('data/crosswind_homotopy.dat','r')
    traj = pickle.load(f)
    f.close()
    for name in [ "y","z",
                  "dy","dz",
                  "w_bn_b_x","w_bn_b_y","w_bn_b_z",
                  "r","dr",'ddr',
                  'ddelta',
                  'aileron','elevator','rudder','flaps',
                  'sin_delta','motor_torque'
                  ]:
        val = traj.lookup(name,timestep=0)
        ocp.bound(name,(val,val),timestep=0,quiet=True)

    dcm0 = rawekite.kiteutils.getDcm(traj, 0, prefix='e')
    dcm1 = rawekite.kiteutils.getDcm(ocp, 0, prefix='e')
    rawekite.kiteutils.matchDcms(ocp, dcm0, dcm1)

    f = open('../pumping_mode/data/crosswind_opt_mechanical_1_loops.dat','r')
    traj = pickle.load(f)
    f.close()
    correspondingNames = {'y':"r_n2b_n_y",
                          'z':"r_n2b_n_z",
                          'dy':"v_bn_n_y",
                          'dz':"v_bn_n_z"}
    obj = 0
    for name in [ 'y',"z",
                  'dy',"dz",
                  "w_bn_b_x","w_bn_b_y","w_bn_b_z",
                  "r","dr",'ddr',
#                  'ddelta',
                  'aileron','elevator','rudder','flaps',
#                  'sin_delta'#,'motor_torque'
                  'e11','e12','e13','e21','e22','e23','e31','e32','e33'
                  ]:
        if name in correspondingNames:
            val = traj.lookup(correspondingNames[name],timestep=-1)
        else:
            val = traj.lookup(name,timestep=-1)
        #ocp.bound(name,(val,val),timestep=-1,quiet=True)
        obj += 1e0*(ocp(name,timestep=-1)-val)**2

    #dcm0 = rawekite.kiteutils.getDcm(traj, -1, prefix='e')
    #dcm1 = rawekite.kiteutils.getDcm(ocp, -1, prefix='e')
    #rawekite.kiteutils.matchDcms(ocp, dcm0, dcm1)

    return (ocp,obj)


if __name__=='__main__':
    from rawe.models.betty_conf import makeConf
    conf = makeConf()
    nk = 120
    dae = carousel_dae.makeDae(conf)
    dae.addP('endTime')

    print "setting up ocp..."
    (ocp,obj) = setupOcp(dae,conf,nk)

    lineRadiusGuess = 4.0

    ocp.setQuadratureDdt('mechanical_energy', 'mechanical_winch_power')
    ocp.setQuadratureDdt('electrical_energy', 'electrical_winch_power')

    # control regularization
    for k in range(ocp.nk):
        dddr = ocp.lookup('dddr',timestep=k)
        daileron = ocp.lookup('daileron',timestep=k)
        delevator = ocp.lookup('delevator',timestep=k)
        drudder = ocp.lookup('drudder',timestep=k)
        dflaps = ocp.lookup('dflaps',timestep=k)
        dmotorTorque = ocp.lookup('dmotor_torque',timestep=k)

        daileronSigma = 0.01
        delevatorSigma = 0.1
        drudderSigma = 0.1
        dflapsSigma = 0.1
        dddrSigma = 0.1
        dmotorTorqueSigma = 0.1

        ailObj = daileron*daileron / (daileronSigma*daileronSigma)
        eleObj = delevator*delevator / (delevatorSigma*delevatorSigma)
        rudObj = drudder*drudder / (drudderSigma*drudderSigma)
        flapsObj = dflaps*dflaps / (dflapsSigma*dflapsSigma)
        winchObj = dddr*dddr / (dddrSigma*dddrSigma)
        motorTorqueObj = dmotorTorque*dmotorTorque / (dmotorTorqueSigma*dmotorTorqueSigma)

        obj += 1e0*(ailObj + eleObj + rudObj + flapsObj + winchObj + motorTorqueObj)/float(ocp.nk)

    ocp.setObjective( obj )

    # initial guesses
    ocp.guess('w0',0)
    ocp.interpolateInitialGuess("data/crosswind_homotopy.dat",force=True,quiet=True,numLoops=5)
    #ocp.interpolateInitialGuess("data/transition.dat",force=True,quiet=True)

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory, showAllPoints=False), 'carousel trajectory')
            ], printBoundViolation=False, printConstraintViolation=False)
#            ], printBoundViolation=True, printConstraintViolation=True)

    # solver
    solverOptions = [("linear_solver","ma27"),
                     ("max_iter",2000),
                     ("expand",True),
                     ("tol",1e-10)]

    print "setting up solver..."
    ocp.setupSolver( solverOpts=solverOptions,
                     callback=callback )

    traj = ocp.solve()
    traj.save("data/transition.dat")

    def printBoundsFeedback():
        # bounds feedback
        xOpt = traj.dvMap.vectorize()
        lbx = ocp.solver.input('lbx')
        ubx = ocp.solver.input('ubx')
        ocp._bounds.printBoundsFeedback(xOpt,lbx,ubx,reportThreshold=0)
    printBoundsFeedback()

    # Plot the results
    def plotResults():
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
