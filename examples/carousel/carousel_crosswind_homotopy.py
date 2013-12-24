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
            #ocp.constrain(ocp.lookup('cL',timestep=k,degIdx=j), '>=', -0.15, tag=('CL > -0.15',(k,j)))
            ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k,degIdx=j), (-10,10), tag=('beta(deg)',(k,j)))
            x = ocp('x', timestep=k,degIdx=j)
            y = ocp('y', timestep=k,degIdx=j)
            z = ocp('z', timestep=k,degIdx=j)
            ocp.constrain(2*C.sqrt(x**2 + y**2), '>=', -z, tag=('azimuth not too high',(k,j)))

def setupOcp(dae,conf,nk,nicp=1,deg=4):
    ocp = rawe.collocation.Coll(dae, nk=nk,nicp=nicp,deg=deg)

    print "setting up collocation..."
    ocp.setupCollocation(ocp.lookup('endTime'))

    # constrain line angle
    for k in range(0,nk):
        ocp.constrain(ocp.lookup('cos_line_angle',timestep=k),'>=',C.cos(65*pi/180), tag=('line angle',k))

    constrainAirspeedAlphaBeta(ocp)
    constrainTetherForce(ocp)
    #realMotorConstraints(ocp)

    # bounds
    ocp.bound('aileron', (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('elevator',(numpy.radians(-10),numpy.radians(10)))
    ocp.bound('rudder',  (numpy.radians(-10),numpy.radians(10)))
    ocp.bound('flaps',  (numpy.radians(0),numpy.radians(0)))
    # can't bound flaps==0 AND have periodic flaps at the same time
    # bounding flaps (-1,1) at timestep 0 doesn't really free them, but satisfies LICQ
    ocp.bound('flaps', (-1,1),timestep=0,quiet=True)
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

    ocp.bound('motor_torque',(-500,500))
    ocp.bound('motor_torque',(0,0),timestep=0)
    #ocp.bound('dmotor_torque',(-1000,1000))
    ocp.bound('dmotor_torque',(0,0))

    ocp.bound('cos_delta',(0,1.5))
    ocp.bound('sin_delta',(-0.4,0.4))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-70,70))

    for w in ['w_bn_b_x',
              'w_bn_b_y',
              'w_bn_b_z']:
        ocp.bound(w,(-4*pi,4*pi))

    ocp.bound('endTime',(3.5,6.0))
#    ocp.bound('endTime',(4.8,4.8))
    ocp.guess('endTime',4.8)
    ocp.bound('w0',(10,10))

    # boundary conditions
    ocp.bound('y',(0,0),timestep=0,quiet=True)
    ocp.bound('sin_delta',(0,0),timestep=0,quiet=True)
    # constrain invariants
    def constrainInvariantErrs():
        rawekite.kiteutils.makeOrthonormal(ocp, ocp.lookup('R_c2b',timestep=0))
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0, tag=('initial c 0',None))
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0, tag=('initial cdot 0',None))
        ocp.constrain(ocp('sin_delta',timestep=0)**2 + ocp('cos_delta',timestep=0)**2,
                      '==', 1, tag=('sin**2 + cos**2 == 1',None))
    constrainInvariantErrs()

    # make it periodic
    for name in [ "x","y","z",
                  "dx","dy","dz",
                  "w_bn_b_x","w_bn_b_y","w_bn_b_z",
                  'ddr',
                  'ddelta',
                  'aileron','elevator','rudder','flaps',
#                  'motor_torque',
                  'sin_delta'
                  ]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1), tag=('periodic '+name,None))

    # periodic attitude
    rawekite.kiteutils.periodicDcm(ocp)

    return ocp


if __name__=='__main__':
    from rawe.models.betty_conf import makeConf
    conf = makeConf()
    conf['runHomotopy'] = True
    nk = 40
    dae = rawe.models.carousel(conf)
    dae.addP('endTime')

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,nk)

    lineRadiusGuess = 80.0
    circleRadiusGuess = 20.0

    # trajectory for homotopy
    homotopyTraj = {'x':[],'y':[],'z':[]}
    # direction =  1: positive about aircraft z
    # direction = -1: negative about aircraft z
    direction = 1
    k = 0
    for nkIdx in range(ocp.nk+1):
        for nicpIdx in range(ocp.nicp):
            if nkIdx == ocp.nk and nicpIdx > 0:
                break
            for degIdx in range(ocp.deg+1):
                if nkIdx == ocp.nk and degIdx > 0:
                    break

                r = circleRadiusGuess
                h = numpy.sqrt(lineRadiusGuess**2 - r**2)
                nTurns = 1

                # path following
                theta = 2*pi*(k+ocp.lagrangePoly.tau_root[degIdx])/float(ocp.nk*ocp.nicp)
                theta -= pi
                theta *= -direction

                thetaDot = nTurns*2*pi/(ocp._guess.lookup('endTime'))
                thetaDot *= -direction
                xyzCircleFrame    = numpy.array([h, r*numpy.sin(theta),          -r*numpy.cos(theta)])
                xyzDotCircleFrame = numpy.array([0, r*numpy.cos(theta)*thetaDot,  r*numpy.sin(theta)*thetaDot])

                phi  = numpy.arcsin(r/lineRadiusGuess) # rotate so it's above ground
                phi += numpy.arcsin((1.3)/lineRadiusGuess)
                phi += 10*pi/180
                R_c2n = numpy.matrix([[  numpy.cos(phi), 0, numpy.sin(phi)],
                                      [               0, 1,              0],
                                      [ -numpy.sin(phi), 0, numpy.cos(phi)]])
                xyz    = numpy.dot(R_c2n, xyzCircleFrame)
                xyzDot = numpy.dot(R_c2n, xyzDotCircleFrame)

                if nicpIdx == 0 and degIdx == 0:
                    homotopyTraj['x'].append(float(xyz[0,0]))
                    homotopyTraj['y'].append(float(xyz[0,1]))
                    homotopyTraj['z'].append(float(xyz[0,2]))

                x = float(xyz[0,0])
                y = float(xyz[0,1])
                z = float(xyz[0,2])

                dx = float(xyzDot[0,0])
                dy = float(xyzDot[0,1])
                dz = float(xyzDot[0,2])

                ocp.guess('x',x,timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('y',y,timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('z',z,timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('dx',dx,timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('dy',dy,timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('dz',dz,timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)

                p0 = numpy.array([x,y,z])
                dp0 = numpy.array([dx,dy,dz])
                e1 = dp0/numpy.linalg.norm(dp0)
                e3 = -p0/lineRadiusGuess
                e2 = numpy.cross(e3,e1)

                ocp.guess('e11',e1[0],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('e12',e1[1],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('e13',e1[2],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)

                ocp.guess('e21',e2[0],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('e22',e2[1],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('e23',e2[2],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)

                ocp.guess('e31',e3[0],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('e32',e3[1],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)
                ocp.guess('e33',e3[2],timestep=nkIdx,nicpIdx=nicpIdx,degIdx=degIdx)

            k += 1
    ocp.guess('w_bn_b_z', direction*2.0*pi/ocp._guess.lookup('endTime'))

    # objective function
    obj = -1e6*ocp.lookup('gamma_homotopy')
    mean_x = numpy.mean(homotopyTraj['x'])
    mean_y = numpy.mean(homotopyTraj['y'])
    mean_z = numpy.mean(homotopyTraj['z'])
    for k in range(ocp.nk+1):
        x = ocp.lookup('x',timestep=k)
        y = ocp.lookup('y',timestep=k)
        z = ocp.lookup('z',timestep=k)
        obj += ((x-mean_x)**2 + (y-mean_y)**2 + (z-mean_z)**2 - circleRadiusGuess**2)**2
        obj += 10*ocp.lookup('sin_delta',timestep=k)**2

#    for k in range(ocp.nk+1):
#        obj += (homotopyTraj['x'][k] - ocp.lookup('x',timestep=k))**2
#        obj += (homotopyTraj['y'][k] - ocp.lookup('y',timestep=k))**2
#        obj += (homotopyTraj['z'][k] - ocp.lookup('z',timestep=k))**2
    ocp.setQuadratureDdt('mechanical_energy', 'mechanical_winch_power')
    ocp.setQuadratureDdt('electrical_energy', 'electrical_winch_power')

    # control regularization
    for k in range(ocp.nk):
        regs = {'dddr':1.0,
                'daileron':numpy.degrees(20.0),
                'delevator':numpy.degrees(20.0),
                'drudder':numpy.degrees(20.0),
                'dflaps':numpy.degrees(20.0),
                'dmotor_torque':5.0}
        for name in regs:
            val = ocp.lookup(name,timestep=k)
            obj += 1e-2*val**2/float(regs[name]**2)/float(ocp.nk)

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
    ocp.guess('w0',10)
    ocp.guess('r',lineRadiusGuess)
    ocp.guess('cos_delta',1)
    ocp.guess('sin_delta',0)

    for name in ['w_bn_b_x','w_bn_b_y','dr','ddr','dddr','aileron','rudder','flaps',
                 'motor_torque','dmotor_torque','ddelta',
                 'elevator','daileron','delevator','drudder','dflaps']:
        ocp.guess(name,0)

    ocp.guess('gamma_homotopy',0)

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory, showAllPoints=True), 'carousel trajectory')
            ], printBoundViolation=True, printConstraintViolation=True)

    # solver
    solverOptions = [("linear_solver","ma97"),
                     ("max_iter",1000),
                     ("expand",True),
                     ("tol",1e-8)]

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

    traj.save("data/carousel_crosswind_homotopy.dat")

    # Plot the results
    def plotResults():
        traj.subplot(['f1_homotopy','f2_homotopy','f3_homotopy'])
        traj.subplot(['t1_homotopy','t2_homotopy','t3_homotopy'])
        traj.subplot(['r_n2b_n_x','r_n2b_n_y','r_n2b_n_z'])
        traj.subplot(['v_bn_n_x','v_bn_n_y','v_bn_n_z'])
        traj.subplot([['aileron','elevator'],['daileron','delevator']],title='control surfaces')
        traj.subplot(['r','dr','ddr'])
        traj.subplot(['wind_at_altitude','dr','v_bn_n_x'])
        traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed')
        traj.subplot([['alpha_deg'],['beta_deg']])
        traj.subplot(['cL','cD','L_over_D'])
        traj.subplot(['mechanical_winch_power', 'tether_tension'])
        traj.subplot(['w_bn_b_x','w_bn_b_y','w_bn_b_z'])
        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        traj.plot(['nu'])
        plt.show()
    plotResults()
