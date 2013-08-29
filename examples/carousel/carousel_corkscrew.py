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
from fourier_fit import TrajFit,FourierFit

import rawe
import rawekite
import carousel_dae
from autogen.tocarouselProto import toProto
from autogen.carousel_pb2 import Trajectory

def constrainTetherForce(ocp):
    for k in range(ocp.nk):
        for j in range(1,ocp.deg+1): #[1]:
            ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=j), '>=', 0, tag=('tether tension positive',k))

def realMotorConstraints(ocp):
    for k in range(nk):
#        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=1),       '<=', 150, tag=('winch torque',k))
#        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=ocp.deg), '<=', 150, tag=('winch torque',k))
        ocp.constrainBnds( ocp.lookup('torque',timestep=k,degIdx=1),
                           (-78,78), tag=('winch torque',k))
        ocp.constrainBnds( ocp.lookup('torque',timestep=k,degIdx=ocp.deg),
                           (-78,78), tag=('winch torque',k))

        ocp.constrain( ocp.lookup('rpm',timestep=k),       '<=', 1500, tag=('rpm',k))
        ocp.constrain( -1500, '<=', ocp.lookup('rpm',timestep=k),       tag=('rpm',k))


def constrainAirspeedAlphaBeta(ocp):
    for k in range(0,ocp.nk):
        for j in range(0,ocp.deg+1):
            ocp.constrainBnds(ocp.lookup('airspeed',timestep=k,degIdx=j),
                              (10,65), tag=('airspeed',(k,j)))
            ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k,degIdx=j), (-4.5,8.5), tag=('alpha(deg)',(k,j)))
            ocp.constrain(ocp.lookup('cL',timestep=k,degIdx=j), '>=', -0.15, tag=('CL > -0.15',(k,j)))
            ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k,degIdx=j), (-6,6), tag=('beta(deg)',(k,j)))
            x = ocp('x', timestep=k,degIdx=j)
            y = ocp('y', timestep=k,degIdx=j)
            z = ocp('z', timestep=k,degIdx=j)
            ocp.constrain(C.sqrt(x**2 + y**2), '>=', -z, tag=('azimuth not too high',(k,j)))

def setupOcp(dae,conf,nk,nicp=1,deg=4):
    ocp = rawe.collocation.Coll(dae, nk=nk,nicp=nicp,deg=deg)

    print "setting up collocation..."
    ocp.setupCollocation(ocp('endTime'))

    # constrain line angle
    for k in range(0,nk):
        ocp.constrain(ocp.lookup('cos_line_angle',timestep=k),'>=',C.cos(55*pi/180), tag=('line angle',k))

    constrainAirspeedAlphaBeta(ocp)
    constrainTetherForce(ocp)
    realMotorConstraints(ocp)

    # bounds
    ocp.bound('aileron', (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('elevator',(numpy.radians(-10),numpy.radians(10)))
    ocp.bound('rudder',  (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('flaps',  (numpy.radians(0),numpy.radians(0)))
    ocp.bound('flaps',  (-1,1), timestep=0, quiet=True) # LICQ, because of IC
    ocp.bound('daileron',  (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('delevator', (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('drudder',   (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('dflaps',    (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('ddelta',(-4*pi,4*pi))

    ocp.bound('x',(-2000,2000))
    ocp.bound('y',(-2000,2000))
    ocp.bound('z',(-2000,1.5))
    ocp.bound('r',(2,300))
    ocp.bound('dr',(-10,10))
    ocp.bound('ddr',(-15,15))
    ocp.bound('dddr',(-50,50))
#    ocp.bound('dr',(-1000,1000))
#    ocp.bound('ddr',(-1500,1500))
#    ocp.bound('dddr',(-500,500))

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

    ocp.bound('endTime',(1,25))
    ocp.bound('w0',(10,10))

    # boundary conditions
    rawekite.kiteutils.makeOrthonormal(ocp, ocp.lookup('R_c2b',timestep=0))
    def getFourierFit(filename,phase):
        # load the fourier fit
        f=open(filename,'r')
        trajFits = pickle.load(f)
        f.close()
        trajFits.setPhase(phase)
        return trajFits

    def get_fourier_dcm(fourier_traj):
        return C.vertcat([C.horzcat([fourier_traj['e11'], fourier_traj['e12'], fourier_traj['e13']]),
                          C.horzcat([fourier_traj['e21'], fourier_traj['e22'], fourier_traj['e23']]),
                          C.horzcat([fourier_traj['e31'], fourier_traj['e32'], fourier_traj['e33']])])

    startup = getFourierFit("data/carousel_homotopy_fourier.dat",ocp('phase0'))

    for name in [ 'x','y','z',
                  'dx','dy','dz',
                  'w_bn_b_x','w_bn_b_y','w_bn_b_z',
                  'r','dr','ddr',
                  'ddelta',
                  'aileron','elevator','rudder','flaps',
                  'cos_delta','sin_delta','motor_torque'
                  ]:
        ocp.constrain(ocp(name,timestep=0),'==',startup[name],tag=('startup '+name,None))

    dcm0 = get_fourier_dcm(startup)
    dcm1 = rawekite.kiteutils.getDcm(ocp, 0, prefix='e')
    rawekite.kiteutils.matchDcms(ocp, dcm0, dcm1, tag='startup')

    crosswind = getFourierFit('../pumping_mode/data/crosswind_opt_mechanical_1_loops_fourier.dat',
                              ocp('phase1'))
    correspondingNames = {'x':'r_n2b_n_x',
                          'y':'r_n2b_n_y',
                          'z':'r_n2b_n_z',
                          'dx':'v_bn_n_x',
                          'dy':'v_bn_n_y',
                          'dz':'v_bn_n_z'}

    obj = 0
    for name in [ 'r','dr','ddr',
                  'x','y','z',
                  'dy','dz',
                  'w_bn_b_x','w_bn_b_y','w_bn_b_z',
#                  'ddelta',
                  'aileron','elevator','rudder','flaps',
#                  'sin_delta','cos_delta',#'motor_torque',
                  'e11','e12','e13','e21','e22','e23','e31','e32','e33'
                  ]:
        if name in correspondingNames:
            name_ = correspondingNames[name]
        else:
            name_ = name
        #ocp.bound(name,(val,val),timestep=-1,quiet=True)
        obj += 1e1*(ocp(name,timestep=-1)-crosswind[name_])**2

    ocp.bound('ddelta',(0,0),timestep=-1,force=True,quiet=True)
    ocp.bound('sin_delta',(0,0),timestep=-1,force=True,quiet=True)
    ocp.bound('cos_delta',(0.5,1.5),timestep=-1,force=True,quiet=True)
#    dcm0 = rawekite.kiteutils.getDcm(traj, -1, prefix='e')
#    dcm1 = rawekite.kiteutils.getDcm(ocp, -1, prefix='e')
#    rawekite.kiteutils.matchDcms(ocp, dcm0, dcm1, tag='crosswind')

    return (ocp,obj)


if __name__=='__main__':
    from rawe.models.betty_conf import makeConf
    conf = makeConf()
    nk = 200
    dae = carousel_dae.makeDae(conf)
    dae.addP('endTime')
    dae.addP('phase0')
    dae.addP('phase1')

    print "setting up ocp..."
    (ocp,obj) = setupOcp(dae,conf,nk)
    ocp.guess('phase0',0)
    ocp.guess('phase1',0)
    ocp.bound('phase0',(-8*pi, 8*pi))
    ocp.bound('phase1',(-8*pi, 8*pi))
#    ocp.bound('phase1',(0,0))

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

        daileronSigma = 0.3
        delevatorSigma = 0.3
        drudderSigma = 0.3
        dflapsSigma = 0.3
        dddrSigma = 5.0
        dmotorTorqueSigma = 10.0

        ailObj = daileron*daileron / (daileronSigma*daileronSigma)
        eleObj = delevator*delevator / (delevatorSigma*delevatorSigma)
        rudObj = drudder*drudder / (drudderSigma*drudderSigma)
        flapsObj = dflaps*dflaps / (dflapsSigma*dflapsSigma)
        winchObj = dddr*dddr / (dddrSigma*dddrSigma)
        motorTorqueObj = dmotorTorque*dmotorTorque / (dmotorTorqueSigma*dmotorTorqueSigma)

        obj += 1e-2*(ailObj + eleObj + rudObj + flapsObj + winchObj + motorTorqueObj)/float(ocp.nk)

    ocp.setObjective( obj )

    # initial guesses
    ocp.interpolateInitialGuess("data/carousel_homotopy.dat",force=True,quiet=True,numLoops=4)
    #ocp.interpolateInitialGuess("data/transition_backup.dat",force=True,quiet=True)

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory, showAllPoints=False), 'carousel trajectory')
            ], printBoundViolation=False, printConstraintViolation=False)
#            ], printBoundViolation=True, printConstraintViolation=True)

    # solver
    solverOptions = [("max_iter",2000),
#                     ("linear_solver","ma86"),
                     ("linear_solver","ma97"),
                     ("expand",True),
#                     ('verbose',True),
                     ("tol",1e-7)]

    print "setting up solver..."
    ocp.setupSolver( solverOpts=solverOptions,
#                     constraintFunOpts=[('verbose',True)],
                     callback=callback)

    import time
    t0 = time.time()
    traj = ocp.solve()
    tTotal = time.time() - t0
    print "total time according to python: " + repr(tTotal)
    traj.save("data/transition.dat")

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
