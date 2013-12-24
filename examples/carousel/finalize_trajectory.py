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
        for j in range(1,ocp.deg+1):
            ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=j), '>=', 0, tag=('tether tension positive',k))

def realMotorConstraints(ocp):
    for k in range(nk-3):
        for j in range(1,ocp.deg+1):
#        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=j),       '<=', 150, tag=('winch torque',(k,j)))
            ocp.constrainBnds( ocp.lookup('torque',timestep=k,degIdx=j),
                               (-78,78), tag=('winch torque',(k,j)))
    for k in range(nk):
        for j in range(0,ocp.deg+1):
            ocp.constrain( ocp.lookup('rpm',timestep=k, degIdx=j),       '<=', 1500, tag=('rpm',k))
            ocp.constrain( -1500, '<=', ocp.lookup('rpm',timestep=k, degIdx=j),       tag=('rpm',k))


def constrainAirspeedAlphaBeta(ocp):
    for k in range(0,ocp.nk):
        for j in range(0,ocp.deg+1):
            ocp.constrainBnds(ocp.lookup('airspeed',timestep=k,degIdx=j),
                              (10,85), tag=('airspeed',(k,j)))
            ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k,degIdx=j), (-4.5,8.5), tag=('alpha(deg)',(k,j)))
            #ocp.constrain(ocp.lookup('cL',timestep=k,degIdx=j), '>=', -0.15, tag=('CL > -0.15',(k,j)))
            ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k,degIdx=j), (-10,10), tag=('beta(deg)',(k,j)))
            #x = ocp('x', timestep=k,degIdx=j)
            #y = ocp('y', timestep=k,degIdx=j)
            #z = ocp('z', timestep=k,degIdx=j)
            #r = C.sqrt(x**2 + y**2)
            #ocp.constrain(-z, '>=', 0.25*r - 2.5, tag=('farther out you go, higher you stay',(k,j)))
            #ocp.constrain(2*C.sqrt(x**2 + y**2), '>=', -z, tag=('azimuth not too high',(k,j)))

def setupOcp(dae,conf,nk,nicp=1,deg=4):
    ocp = rawe.collocation.Coll(dae, nk=nk,nicp=nicp,deg=deg)

    print "setting up collocation..."
    ocp.setupCollocation(ocp('endTime'))

    # constrain line angle
    for k in range(0,nk):
        ocp.constrain(ocp.lookup('cos_line_angle',timestep=k),'>=',C.cos(65*pi/180), tag=('line angle',k))

    constrainAirspeedAlphaBeta(ocp)
    constrainTetherForce(ocp)
#    realMotorConstraints(ocp)

    # bounds
    ocp.bound('aileron', (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('elevator',(numpy.radians(-10),numpy.radians(10)))
    ocp.bound('rudder',  (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('flaps',  (numpy.radians(0),numpy.radians(0)))
    ocp.bound('flaps',  (-1,1), timestep=0, quiet=True) # LICQ, because of IC
    ocp.bound('flaps',  (-1,1), timestep=-1, quiet=True) # LICQ, because of FC
    ocp.bound('daileron',  (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('delevator', (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('drudder',   (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('dflaps',    (numpy.radians(-40), numpy.radians(40)))
    ocp.bound('ddelta',(-2*pi, 2*pi))

    ocp.bound('x',(-2000,2000))
    ocp.bound('y',(-2000,2000))
    ocp.bound('z',(-2000,0))
    ocp.bound('r',(2,300))
    ocp.bound('dr',(-100,100))
    ocp.bound('ddr',(-150,150))
    ocp.bound('dddr',(-200,200))
#    ocp.bound('dr',(-1000,1000))
#    ocp.bound('ddr',(-1500,1500))
#    ocp.bound('dddr',(-500,500))

    ocp.bound('motor_torque',(-300,300))
    ocp.bound('dmotor_torque',(-10000,10000))

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

    ocp.bound('endTime',(1,24))
    ocp.bound('w0',(10,10))

    # boundary conditions
    # constrain invariants
    def constrainInvariantErrs():
        rawekite.kiteutils.makeOrthonormal(ocp, ocp.lookup('R_c2b',timestep=0))
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0, tag=('initial c 0',None))
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0, tag=('initial cdot 0',None))
        ocp.constrain(ocp('sin_delta',timestep=0)**2 + ocp('cos_delta',timestep=0)**2,
                      '==', 1, tag=('sin**2 + cos**2 == 1',None))
    constrainInvariantErrs()

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
                  'ddr',
                  'ddelta',
                  'aileron','elevator','rudder','flaps',
                  'motor_torque'
                  ]:
        ocp.constrain(ocp(name,timestep=0),'==',startup[name],tag=('initial '+name,None))
    ocp.constrain(C.arctan2(ocp('sin_delta',timestep=0), ocp('cos_delta',timestep=0)),
                  '==',
                  C.arctan2(startup['sin_delta'], startup['cos_delta']),
                  tag=('startup delta matching',None))

    dcm0 = get_fourier_dcm(startup)
    dcm1 = rawekite.kiteutils.getDcm(ocp, 0, prefix='e')
    rawekite.kiteutils.matchDcms(ocp, dcm0, dcm1, tag='startup')

    crosswind = getFourierFit('../pumping_mode/data/crosswind_opt_mechanical_6_loops_fourier.dat',
                              ocp('phase1'))
    correspondingNames = {'x':'r_n2b_n_x',
                          'y':'r_n2b_n_y',
                          'z':'r_n2b_n_z',
                          'dx':'v_bn_n_x',
                          'dy':'v_bn_n_y',
                          'dz':'v_bn_n_z'}

    for name in [ 'x','y','z',
                  'dx','dy','dz',
                  'ddr',
                  'w_bn_b_x','w_bn_b_y','w_bn_b_z',
#                  'aileron','elevator','rudder','flaps',
                  ]:
        if name in correspondingNames:
            name_ = correspondingNames[name]
        else:
            name_ = name
        ocp.constrain(ocp(name,timestep=-1), '==', crosswind[name_], tag=('terminal '+name_,None))

    ocp.bound('ddelta',(0,0),timestep=-1,force=True,quiet=True)
    ocp.bound('sin_delta',(0,0),timestep=-1,force=True,quiet=True)
    ocp.bound('cos_delta',(0.5,1.5),timestep=-1,force=True,quiet=True)
    dcm0 = get_fourier_dcm(crosswind)
    dcm1 = rawekite.kiteutils.getDcm(ocp, -1, prefix='e')
    rawekite.kiteutils.matchDcms(ocp, dcm0, dcm1, tag='crosswind')

    return ocp


if __name__=='__main__':
    from rawe.models.betty_conf import makeConf
    conf = makeConf()
    nk = 200
    dae = carousel_dae.makeDae(conf)
    dae.addP('endTime')
    dae.addP('phase0')
    dae.addP('phase1')

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,nk)
    ocp.bound('phase0',(-8*pi, 8*pi))
    #ocp.bound('phase1',(0.05*2*pi, 2*pi))
    ocp.bound('phase1',(1.8, 2*pi))
    #ocp.bound('phase1',(numpy.radians(45), 8*pi))

    ocp.setQuadratureDdt('mechanical_energy', 'mechanical_winch_power')
    ocp.setQuadratureDdt('electrical_energy', 'electrical_winch_power')

    # control regularization
    reg_obj = 0
    for k in range(ocp.nk):
        regs = {'dddr':1.0,
                'daileron':numpy.degrees(20.0),
                'delevator':numpy.degrees(20.0),
                'drudder':numpy.degrees(20.0),
                'dflaps':numpy.degrees(20.0),
                'dmotor_torque':5.0}
        for name in regs:
            val = ocp.lookup(name,timestep=k)
            reg_obj += val**2/float(regs[name]**2)/float(ocp.nk)

    # state regularization
    for k in range(ocp.nk):
        for j in range(ocp.deg+1):
            regs = {'w_bn_b_x':1.0,
                    'w_bn_b_y':1.0,
                    'w_bn_b_z':1.0,
                    'dr':5.0,
                    'ddr':10.0,
                    'motor_torque':20.0,
                    'aileron':numpy.degrees(10.0),
                    'elevator':numpy.degrees(10.0),
                    'rudder':numpy.degrees(10.0),
                    'flaps':numpy.degrees(10.0),
                    'beta_deg':4.0,
                    'alpha_deg':25.0}
            for name in regs:
                val = ocp.lookup(name,timestep=k,degIdx=j)
                reg_obj += val**2/float(regs[name]**2)/float(ocp.nk*(ocp.deg+1))

    ocp.setObjective( reg_obj )

    # initial guesses
    #ocp.interpolateInitialGuess("data/transition.dat",force=True,quiet=True)
    ocp.interpolateInitialGuess("data/transition_backup.dat",force=True,quiet=True)
    #ocp.interpolateInitialGuess("data/final_transition_backup.dat",force=True,quiet=True)

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory, showAllPoints=False), 'carousel trajectory')
            #], printBoundViolation=False, printConstraintViolation=False)
        ], printBoundViolation=True, printConstraintViolation=True)

    # solver
    solverOptions = [("max_iter",2000),
#                     ("linear_solver","ma86"),
                     ("linear_solver","ma97"),
                     ("expand",True),
#                     ('verbose',True),
                     ("tol",1e-4)]

    print "setting up solver..."
    ocp.setupSolver( solverOpts=solverOptions,
#                     constraintFunOpts=[('verbose',True)],
                     callback=callback)

    import time
    t0 = time.time()
    traj = ocp.solve()
    tTotal = time.time() - t0
    print "total time according to python: " + repr(tTotal)
    traj.save("data/final_transition.dat")

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
#    plotResults()
