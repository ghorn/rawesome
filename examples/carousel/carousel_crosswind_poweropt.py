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

import casadi as C
import matplotlib.pyplot as plt
import pickle
import numpy
from numpy import pi

import rawe
import rawekite
import carousel_dae
from autogen.tocarouselProto import toProto
from autogen.carousel_pb2 import Trajectory

numLoops=1
powerType = 'mechanical'
#powerType = 'electrical'

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
                              (10,75), tag=('airspeed',(k,j)))
            ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k,degIdx=j), (-4.5,8.5), tag=('alpha(deg)',(k,j)))
            #ocp.constrain(ocp.lookup('cL',timestep=k,degIdx=j), '>=', -0.15, tag=('CL > -0.15',(k,j)))
            ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k,degIdx=j), (-10,10), tag=('beta(deg)',(k,j)))

def setupOcp(dae,conf,nk,nicp=1,deg=4):
    ocp = rawe.collocation.Coll(dae, nk=nk, nicp=nicp, deg=deg)

    print "setting up collocation..."
    ocp.setupCollocation(ocp.lookup('endTime'))

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
    ocp.bound('daileron',  (numpy.radians(-80), numpy.radians(80)))
    ocp.bound('delevator', (numpy.radians(-80), numpy.radians(80)))
    ocp.bound('drudder',   (numpy.radians(-80), numpy.radians(80)))
    ocp.bound('dflaps',    (numpy.radians(-80), numpy.radians(80)))
    ocp.bound('ddelta',(-2*pi, 2*pi))

    ocp.bound('x',(0,2000))
    ocp.bound('y',(-2000,2000))
    ocp.bound('z',(-2000,2))
    ocp.bound('r',(2,300))
    ocp.bound('dr', (-100,100))
        
    ocp.bound('ddr',(-1500,1500))
    ocp.bound('dddr',(-50000,50000))
#    ocp.bound('dr',(-1000,1000))
#    ocp.bound('ddr',(-1500,1500))
#    ocp.bound('dddr',(-500,500))

    ocp.bound('motor_torque',(-300,300))
    #ocp.bound('motor_torque',(0,0))
    #ocp.bound('motor_torque',(0,0),timestep=0)
    ocp.bound('dmotor_torque',(-100000,100000))
    #ocp.bound('dmotor_torque',(0,0))

    ocp.bound('cos_delta',(0,1.5))
    ocp.bound('sin_delta',(-1.5,1.5))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-1000,1000))

    for w in ['w_bn_b_x',
              'w_bn_b_y',
              'w_bn_b_z']:
        ocp.bound(w,(-6*pi,6*pi))

    # boundary conditions
    # constrain invariants
    def constrainInvariantErrs():
        rawekite.kiteutils.makeOrthonormal(ocp, ocp.lookup('R_c2b',timestep=0))
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0, tag=('initial c 0',None))
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0, tag=('initial cdot 0',None))
        ocp.constrain(ocp('sin_delta',timestep=0)**2 + ocp('cos_delta',timestep=0)**2,
                      '==', 1, tag=('sin**2 + cos**2 == 1',None))
    constrainInvariantErrs()

    # make it periodic
    for name in [ 'x','y','z',
                  'dx','dy','dz',
                  'w_bn_b_x','w_bn_b_y','w_bn_b_z',
                  'ddr',
                  'ddelta',
                  'aileron','elevator','rudder','flaps',
                  'motor_torque',
                  'sin_delta'
                  ]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1), tag=('periodic diff state \"'+name+'"',None))

    # periodic attitude
    rawekite.kiteutils.periodicDcm(ocp)

#    ocp.bound('endTime',(0.5,12))
    ocp.bound('endTime',(1.0,numLoops*7.5))
    ocp.bound('w0',(10,10))

    # boundary conditions
    #ocp.bound('y',(0,0),timestep=0,quiet=True)
    ocp.bound('sin_delta',(0,0),timestep=0,quiet=True)

    # objective function
    # control regularization
    reg_obj = 0
    for k in range(ocp.nk):
        regs = {'dddr':1000.0,
                'daileron':numpy.degrees(20.0),
                'delevator':numpy.degrees(20.0),
                'drudder':numpy.degrees(20.0),
                'dflaps':numpy.degrees(20.0),
                'dmotor_torque':3000.0}
        for name in regs:
            val = ocp.lookup(name,timestep=k)
            reg_obj += val**2/float(regs[name]**2)/float(ocp.nk)

    # state regularization
    for k in range(ocp.nk):
        for j in range(ocp.deg+1):
            regs = {'w_bn_b_x':pi,
                    'w_bn_b_y':pi,
                    'w_bn_b_z':pi,
                    'dr':100.0,
#                    'ddr':10.0,
                    'motor_torque':300.0,
                    'aileron':numpy.degrees(10.0),
                    'elevator':numpy.degrees(10.0),
                    'rudder':numpy.degrees(10.0),
                    'flaps':numpy.degrees(10.0),
                    'beta_deg':5.0,
                    'alpha_deg':25.0}
            for name in regs:
                val = ocp.lookup(name,timestep=k,degIdx=j)
                reg_obj += val**2/float(regs[name]**2)/float(ocp.nk*(ocp.deg+1))

    ocp.setQuadratureDdt('electrical_winch_energy', 'electrical_winch_power')
    ocp.setQuadratureDdt('mechanical_winch_energy', 'mechanical_winch_power')
    ocp.setQuadratureDdt('mechanical_rot_motor_energy', 'mechanical_rot_motor_power')
    #ocp.setQuadratureDdt('mechanical_rot_motor_energy', 'mechanical_rot_motor_power_slack')
    #for k in range(ocp.nk):
    #    for j in range(1,ocp.deg+1):
    #        for i in range(ocp.nicp):
    #            slack = ocp('mechanical_rot_motor_power_slack',timestep=k,degIdx=j,nicpIdx=i)
    #            val   = ocp('mechanical_rot_motor_power',timestep=k,degIdx=j,nicpIdx=i)
    #            ocp.constrain(val, '<=', slack, tag=('rot motor slack',(k,j,i)))
    #            ocp.constrain(0,   '<=', slack, tag=('rot motor slack 0',(k,j,i)))

    energy = ocp.lookup(powerType+'_winch_energy',timestep=-1)
    energy += ocp.lookup('mechanical_rot_motor_energy',timestep=-1)
    ocp.setObjective( 1e2*reg_obj + energy/ocp.lookup('endTime') )

    return ocp


if __name__=='__main__':
    from rawe.models.betty_conf import makeConf
    conf = makeConf()
#    nk = 30*numLoops
    nk = 60
    dae = carousel_dae.makeDae(conf)
    dae.addP('endTime')
    #dae.addZ('mechanical_rot_motor_power_slack')

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,nk)

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory, showAllPoints=False), 'carousel trajectory')
        ], printBoundViolation=True, printConstraintViolation=True)

    # solver
    solverOptions = [("max_iter",2000),
                     ("linear_solver","ma97"),
                     ("expand",True),
                     ("tol",1e-8)]
    print "setting up solver..."
    ocp.setupSolver( solverOpts=solverOptions,
                     callback=callback )
    ocp.interpolateInitialGuess("data/carousel_crosswind_homotopy.dat",force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess("data/carousel_crosswind_opt_mechanical_1_loops.dat",force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess("data/crosswind_opt_mechanical_6_loops-backup.dat",force=True,quiet=True)
#    ocp.interpolateInitialGuess("../pumping_mode/data/crosswind_opt_mechanical_6_loops-backup.dat",force=True,quiet=True)
#    ocp.interpolateInitialGuess("data/crosswind_opt_electrical_1_loops.dat",force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess('data/crosswind_opt_'+powerType+'_1_loops.dat',
#                                force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess("data/crosswind_opt.dat",force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess("data/crosswind_opt_electrical_2_loops.dat",force=True,quiet=True,numLoops=1)

    import time
    t0 = time.time()
    traj = ocp.solve()
    tTotal = time.time() - t0
    print "total time according to python: " + repr(tTotal)

    print "num loops: "+str(numLoops)
    print "optimizing "+powerType
    print "optimal electrical winch power: "+str(traj.lookup('electrical_winch_energy',-1)/traj.lookup('endTime'))
    mech_e = traj.lookup('mechanical_winch_energy',-1)/traj.lookup('endTime')
#    mech_rot_e = traj.lookup('mechanical_rot_motor_energy',-1)/traj.lookup('endTime')
    print "optimal mechanical winch power: "+str(mech_e)
    print "optimal mechanical rot motor power: "+str(mech_rot_e)
    print "endTime: "+str(traj.lookup('endTime'))
#    print "total power: "+str(mech_e + mech_rot_e)

#    traj.saveMat('data/carousel_crosswind_opt_'+powerType+'_'+str(numLoops)+'_loops.mat',
#                 dataname='carousel_crosswind_opt_'+powerType+'_'+str(numLoops)+'_loops')
    traj.save('data/carousel_crosswind_opt_'+powerType+'_'+str(numLoops)+'_loops.dat')

    # Plot the results
    def plotResults():
#        traj.subplot(['aero_fx','aero_fy','aero_fz'])
#        traj.subplot(['aero_mx','aero_my','aero_mz'])
#        traj.subplot(['r_n2b_n_x','r_n2b_n_y','r_n2b_n_z'])
#        traj.subplot(['v_bn_n_x','v_bn_n_y','v_bn_n_z'])
        traj.subplot(['aileron','elevator','rudder','flaps'],title='control surfaces')
        traj.subplot([['dddr'],['daileron','delevator','drudder','dflaps']],title='control surfaces')
#        traj.subplot(['wind_at_altitude','dr'],title='')
#        traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed',title='airspeed')
        traj.subplot([['alpha_deg'],['beta_deg']])
        traj.subplot([['cL'],['cD','cD_tether'],['L_over_D','L_over_D_with_tether']],title='')
#        traj.subplot([['winch_power'], ['tether_tension'],['accel_g','accel_without_gravity_g']])
#        traj.subplot([['rpm'],['dr']])
        traj.subplot([['tether_tension'],['torque']])
#        traj.plot(['mechanical_winch_power', 'electrical_winch_power'])
#        traj.plot('r')
#        traj.subplot([['ddx','ddy','ddz'],['accel','accel without gravity']])
#        traj.plot(["loyds_limit","loyds_limit_exact","neg_winch_power"])
#        traj.plot(["loyd's limit","-(winch power)"],title='')
        traj.subplot(['daileronCost','delevatorCost','dddrCost'])
        traj.subplot(['r','dr','ddr','dddr'])
#        traj.subplot(['w_bn_b_x','w_bn_b_y','w_bn_b_z'])
#        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
#        traj.plot('line_angle_deg')
#        traj.plot('quadrature_energy')
#        traj.subplot(['energy','quadrature_energy'])
#        traj.plot(['energy','quadrature_energy'])
#        traj.plot('nu')

        plt.show()
#    plotResults()
#    traj.subplot(['r','dr'])
#    traj.subplot(['rpm','torque'])
#    traj.subplot(['airspeed'])
#    plt.show()
