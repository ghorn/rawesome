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
from autogen.topumpingProto import toProto
from autogen.pumping_pb2 import Trajectory

numLoops=1
powerType = 'mechanical'
#powerType = 'electrical'

def constrainInvariants(ocp):
    R_n2b = ocp.lookup('R_n2b',timestep=0)
    rawekite.kiteutils.makeOrthonormal(ocp, R_n2b)
    ocp.constrain(ocp.lookup('c',timestep=0), '==', 0, tag=('c(0)==0',None))
    ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0, tag=('cdot(0)==0',None))

    # constrain line angle
    for k in range(0,nk):
        for j in range(0,ocp.deg+1):
            ocp.constrain(ocp.lookup('cos_line_angle',timestep=k,degIdx=j),'>=',C.cos(65*pi/180), tag=('line angle',k))

def constrainAirspeedAlphaBeta(ocp):
    for k in range(0,nk):
        ocp.constrain(ocp.lookup('airspeed',timestep=k), '>=', 10, tag=('airspeed',k))
        for j in range(0,ocp.deg+1):
            ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k,degIdx=j), (-4.5,8.5), tag=('alpha(deg)',k))

        ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k), (-9,9), tag=('beta(deg)',k))

def constrainTetherForce(ocp):
    for k in range(ocp.nk):
#        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=1), '>=', 0, tag=('tether tension positive',k))
#        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=ocp.deg), '>=', 0, tag=('tether tension positive',k))
        for j in range(1,ocp.deg+1):
            ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=j), '>=', 0, tag=('tether tension positive',k))

def realMotorConstraints(ocp):
    for k in range(nk):
#        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=1),       '<=', 150, tag=('motor torque',k))
#        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=ocp.deg), '<=', 150, tag=('motor torque',k))
        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=1),       '<=', 78, tag=('motor torque',k))
        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=ocp.deg), '<=', 78, tag=('motor torque',k))

        ocp.constrain( ocp.lookup('rpm',timestep=k),       '<=', 1500, tag=('rpm',k))
        ocp.constrain( -1500, '<=', ocp.lookup('rpm',timestep=k),       tag=('rpm',k))


def setupOcp(dae,conf,nk,nicp,deg,collPoly):
    def addCosts():
        dddr = dae['dddr']
        daileron = dae['daileron']
        delevator = dae['delevator']
        drudder = dae['drudder']
        dflaps = dae['dflaps']

        daileronSigma = 0.001
        delevatorSigma = 0.8
        dddrSigma = 10.0
        drudderSigma = 0.1
        dflapsSigma = 0.1

        fudgeFactor = 1e-1
        nkf = float(nk)
        dae['daileronCost'] =  fudgeFactor*daileron*daileron / (daileronSigma*daileronSigma*nkf)
        dae['delevatorCost'] = fudgeFactor*delevator*delevator / (delevatorSigma*delevatorSigma*nkf)
        dae['drudderCost'] =   fudgeFactor*drudder*drudder / (drudderSigma*drudderSigma*nkf)
        dae['dflapsCost'] =    fudgeFactor*dflaps*dflaps / (dflapsSigma*dflapsSigma*nkf)
        dae['dddrCost'] =      fudgeFactor*dddr*dddr / (dddrSigma*dddrSigma*nkf)

    addCosts()

    ocp = rawe.collocation.Coll(dae, nk=nk, nicp=nicp, deg=deg, collPoly=collPoly)

    ocp.setupCollocation(ocp.lookup('endTime'))

    constrainInvariants(ocp)
    constrainAirspeedAlphaBeta(ocp)
    constrainTetherForce(ocp)
    realMotorConstraints(ocp)

    # make it periodic
    for name in [ "r_n2b_n_y","r_n2b_n_z",
                  "v_bn_n_y","v_bn_n_z",
                  "w_bn_b_x","w_bn_b_y","w_bn_b_z",
                  "r","dr","ddr",
                  'aileron','elevator','rudder','flaps'
                  ]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1), tag=('periodic diff state \"'+name+'"',None))

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
    ocp.bound('daileron',(-2.0,2.0))
    ocp.bound('delevator',(-2.0,2.0))
    ocp.bound('drudder',(-2.0,2.0))
    ocp.bound('dflaps',(-2.0,2.0))

    ocp.bound('r_n2b_n_x',(-2000,2000))
    ocp.bound('r_n2b_n_y',(-2000,2000))
    if 'minAltitude' in conf:
        ocp.bound('r_n2b_n_z',(-2000, -conf['minAltitude']))
    else:
        ocp.bound('r_n2b_n_z',(-2000, -0.05))
    ocp.bound('r',(1,500))
    ocp.bound('dr',(-30,30))
    ocp.bound('ddr',(-500,500))
    ocp.bound('dddr',(-50000,50000))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['v_bn_n_x','v_bn_n_y','v_bn_n_z']:
        ocp.bound(d,(-200,200))

    for w in ['w_bn_b_x',
              'w_bn_b_y',
              'w_bn_b_z']:
        ocp.bound(w,(-6*pi,6*pi))

#    ocp.bound('endTime',(0.5,12))
    ocp.bound('endTime',(0.5,numLoops*7.5))
    ocp.bound('w0',(10,10))

    # boundary conditions
    ocp.bound('r_n2b_n_y',(0,0),timestep=0,quiet=True)

    # guesses
    ocp.guess('endTime',5.4)
    ocp.guess('w0',10)

    # objective function
    obj = 0
    for k in range(nk):
        # control regularization
        obj += ocp.lookup('daileronCost',timestep=k)
        obj += ocp.lookup('delevatorCost',timestep=k)
        obj += ocp.lookup('drudderCost',timestep=k)
        obj += ocp.lookup('dflapsCost',timestep=k)
        obj += ocp.lookup('dddrCost',timestep=k)

    ocp.setQuadratureDdt('mechanical_energy', 'mechanical_winch_power')
    ocp.setQuadratureDdt('electrical_energy', 'electrical_winch_power')

    ocp.setObjective( obj + ocp.lookup(powerType+'_energy',timestep=-1)/ocp.lookup('endTime') )

    return ocp


if __name__=='__main__':
    print "reading config..."
#    from carousel_conf import conf
    #from highwind_carousel_conf import conf
    from rawe.models.betty_conf import makeConf

    nk = 100*numLoops
#    nk = 70

    print "creating model..."
    conf = makeConf()
    dae = rawe.models.crosswind(conf)
    dae.addP('endTime')
    conf['minAltitude'] = 0

    print "setting up ocp..."
    nicp = 1
    deg = 4
    collPoly='RADAU'
    #collPoly='LEGENDRE'
    ocp = setupOcp(dae,conf,nk,nicp,deg,collPoly)

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(
        ocp, callbacks=[
            (rawe.telemetry.trajectoryCallback(toProto, Trajectory, showAllPoints=False), 'pumping trajectory')
        ])

    # solver
    ipoptOptions = [("linear_solver","ma97"),
                    ("expand",True),
                    ("max_iter",2000),
                    ("tol",1e-8)]
    worhpOptions = [("Max_Iter",5000),
                    ("expand",True),
                    #("MaxIter",5000),
                    ("Timeout", 1e6),
                    ("UserHM", True),
                    ("ScaleConIter",True),
                    ("ScaledFD",True),
                    ("ScaledKKT",True),
                    ("ScaledObj",True),
                    ("ScaledQP",True)
                    ]
    print "setting up solver..."
    solverOptions = ipoptOptions
#    solverOptions = worhpOptions
    ocp.setupSolver( solverOpts=solverOptions,
                     callback=callback )

    ocp.interpolateInitialGuess("data/crosswind_homotopy.dat",force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess("data/crosswind_opt_electrical_1_loops.dat",force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess('data/crosswind_opt_'+powerType+'_1_loops.dat',
#                                force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess("data/crosswind_opt.dat",force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess("data/crosswind_opt_electrical_2_loops.dat",force=True,quiet=True,numLoops=1)

    traj = ocp.solve()

    print "num loops: "+str(numLoops)
    print "optimizing "+powerType
    print "optimal mechanical power: "+str(traj.lookup('mechanical_energy',-1)/traj.lookup('endTime'))
    print "optimal electrical power: "+str(traj.lookup('electrical_energy',-1)/traj.lookup('endTime'))
    print "endTime: "+str(traj.lookup('endTime'))

    traj.saveMat('data/crosswind_opt_'+powerType+'_'+str(numLoops)+'_loops.mat',
                 dataname='crosswind_opt_'+powerType+'_'+str(numLoops)+'_loops')
    traj.save('data/crosswind_opt_'+powerType+'_'+str(numLoops)+'_loops.dat')

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
    traj.subplot(['r','dr','ddr','dddr'])
    plt.show()


    def plotPaper():
        plt.figure()
        plt.subplot(221)
        traj._plot('cL','',showLegend=False)
        plt.legend(['$C_L$'])
        plt.xlabel('')
        plt.ylim([0.1,1.1])
        plt.xlim([0,traj.tgrid[-1,0,0]])

        plt.subplot(223)
        traj._plot('L_over_D_with_tether','',showLegend=False)
        plt.legend(['$L/D$'])
        plt.ylim([2,15])
        plt.xlim([0,traj.tgrid[-1,0,0]])

        plt.subplot(222)
        traj._plot('wind_at_altitude','',showLegend=False)
        plt.xlabel('')
        plt.ylabel('[m/s]')
        plt.ylim([5.5,10])
        plt.legend(['wind at altitude'])
        plt.xlim([0,traj.tgrid[-1,0,0]])

        plt.subplot(224)
        traj._plot('dr','',showLegend=False)
        plt.ylabel('[m/s]')
        plt.ylim([-15,12])
        plt.legend(['$\dot{l}$'])
        plt.xlim([0,traj.tgrid[-1,0,0]])

#        traj.subplot(['cL','L/D','dr'],title='')
#        traj.plot(["loyd's limit","loyd's limit (exact)","-(winch power)"])
#        traj.plot(["loyd's limit","-(winch power)"],title='')
        plt.figure()
        traj._plot("loyds_limit",'',showLegend=False,style=':')
        traj._plot("neg_winch_power",'',showLegend=False)
        plt.legend(["Loyd's limit","winch power"])
        plt.ylabel('power [W]')
        plt.ylim([-600,1100])
        plt.grid()
        plt.xlim([0,traj.tgrid[-1,0,0]])


        plt.show()
    #plotPaper()
