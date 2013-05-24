import casadi as C
import matplotlib.pyplot as plt
import pickle
import numpy
from numpy import pi

import rawe
import rawekite

numLoops=1
#powerType = 'electrical'
powerType = 'mechanical'

def setupOcp(dae,conf,nk,nicp,deg,collPoly):
    def addCosts():
        dddr = dae['dddr']
        daileron = dae['daileron']
        delevator = dae['delevator']

        daileronSigma = 0.001
        delevatorSigma = 0.8
        dddrSigma = 20.0

        fudgeFactor = 1e-1
        nkf = float(nk)
        dae['daileronCost'] =  fudgeFactor*daileron*daileron / (daileronSigma*daileronSigma*nkf)
        dae['delevatorCost'] = fudgeFactor*delevator*delevator / (delevatorSigma*delevatorSigma*nkf)
        dae['dddrCost'] =      fudgeFactor*dddr*dddr / (dddrSigma*dddrSigma*nkf)

    addCosts()

    ocp = rawe.collocation.Coll(dae, nk=nk, nicp=nicp, deg=deg, collPoly=collPoly)

    print "setting up collocation..."
    ocp.setupCollocation(ocp.lookup('endTime'))

    print "moar setting up ocp..."

    # constrain invariants
    def constrainInvariantErrs():
        dcm = ocp.lookup('dcm',timestep=0)
        rawekite.kiteutils.makeOrthonormal(ocp, dcm)
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0, tag=('c(0)==0',None))
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0, tag=('cdot(0)==0',None))
    constrainInvariantErrs()

    # constrain line angle
    for k in range(0,nk):
        for j in range(0,ocp.deg+1):
            ocp.constrain(ocp.lookup('cos_line_angle',timestep=k,degIdx=j),'>=',C.cos(65*pi/180), tag=('line angle',k))

    # constrain airspeed
    def constrainAirspeedAlphaBeta():
        for k in range(0,nk):
            ocp.constrain(ocp.lookup('airspeed',timestep=k), '>=', 10, tag=('airspeed',k))
            for j in range(0,ocp.deg+1):
                ocp.constrainBnds(ocp.lookup('alpha_deg',timestep=k,degIdx=j), (-8.5,9.5), tag=('alpha(deg)',k))

            ocp.constrainBnds(ocp.lookup('beta_deg', timestep=k), (-10,10), tag=('beta(deg)',k))
    constrainAirspeedAlphaBeta()
    #def constrainCl():
    #    for k in range(0,nk):
    #        ocp.constrain(ocp.lookup('cL',timestep=k), '<=', 2.0, tag=('cL',k))
    #constrainCl()

    # constrain tether force
    for k in range(nk):
#        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=1), '>=', 0, tag=('tether tension positive',k))
#        ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=ocp.deg), '>=', 0, tag=('tether tension positive',k))
        for j in range(1,ocp.deg+1):
            ocp.constrain( ocp.lookup('tether_tension',timestep=k,degIdx=j), '>=', 0, tag=('tether tension positive',k))

    # real motor constraints
    for k in range(nk):
#        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=1),       '<=', 150, tag=('motor torque',k))
#        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=ocp.deg), '<=', 150, tag=('motor torque',k))
        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=1),       '<=', 78, tag=('motor torque',k))
        ocp.constrain( ocp.lookup('torque',timestep=k,degIdx=ocp.deg), '<=', 78, tag=('motor torque',k))

        ocp.constrain( ocp.lookup('rpm',timestep=k),       '<=', 1500, tag=('rpm',k))
        ocp.constrain( -1500, '<=', ocp.lookup('rpm',timestep=k),       tag=('rpm',k))

    # make it periodic
    for name in [ "y","z",
                  "dy","dz",
                  "w_bn_b_x","w_bn_b_y","w_bn_b_z",
                  "r","dr","ddr",
                  'aileron','elevator'
                  ]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1), tag=('periodic diff state \"'+name+'"',None))

    # periodic attitude
    rawekite.kiteutils.periodicDcm(ocp)

    # bounds
    ocp.bound('aileron',(-0.04,0.04))
    ocp.bound('elevator',(-0.1,0.5))
    ocp.bound('daileron',(-2.0,2.0))
    ocp.bound('delevator',(-2.0,2.0))

    ocp.bound('x',(-2000,2000))
    ocp.bound('y',(-2000,2000))
    if 'minAltitude' in conf:
        ocp.bound('z',(conf['minAltitude'],2000))
    else:
        ocp.bound('z',(0.01,2000))
    ocp.bound('r',(1,500))
    ocp.bound('dr',(-30,30))
    ocp.bound('ddr',(-500,500))
    ocp.bound('dddr',(-50000,50000))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-200,200))

    for w in ['w_bn_b_x',
              'w_bn_b_y',
              'w_bn_b_z']:
        ocp.bound(w,(-6*pi,6*pi))

#    ocp.bound('endTime',(0.5,12))
    ocp.bound('endTime',(0.5,numLoops*7.5))
    ocp.bound('w0',(10,10))

    # boundary conditions
    ocp.bound('y',(0,0),timestep=0,quiet=True)

    # guesses
    ocp.guess('endTime',5.4)
    ocp.guess('w0',10)

    # objective function
    obj = 0
    for k in range(nk):
        # control regularization
        obj += ocp.lookup('daileronCost',timestep=k)
        obj += ocp.lookup('delevatorCost',timestep=k)
        obj += ocp.lookup('dddrCost',timestep=k)

    ocp.setQuadratureDdt('mechanical_energy', 'mechanical_winch_power')
    ocp.setQuadratureDdt('electrical_energy', 'electrical_winch_power')

    ocp.setObjective( obj + ocp.lookup(powerType+'_energy',timestep=-1)/ocp.lookup('endTime') )

    return ocp


if __name__=='__main__':
    print "reading config..."
#    from carousel_conf import conf
    #from highwind_carousel_conf import conf
    from betty_conf import conf

    nk = 60*numLoops
#    nk = 70

    print "creating model..."
    dae = rawe.models.crosswind(conf)
    dae.addP('endTime')

    print "setting up ocp..."
    nicp = 1
    deg = 4
    collPoly='RADAU'
    #collPoly='LEGENDRE'
    ocp = setupOcp(dae,conf,nk,nicp,deg,collPoly)

    # spawn telemetry thread
    callback = rawe.telemetry.startTelemetry(ocp, conf, callbacks=[(rawekite.kiteTelemetry.showAllPoints,'multi-carousel')])
#    callback = rawe.telemetry.startTelemetry(ocp, conf, callbacks=[(rawekite.kiteTelemetry.showAllPoints,'multi-carousel')], printBoundViolation=True, printConstraintViolation=True)

    # solver
    ipoptOptions = [("linear_solver","ma27"),
                    ("max_iter",1000),
                    ("tol",1e-12)]
    worhpOptions = [("Max_Iter",5000),
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
#    ocp.interpolateInitialGuess("data/crosswind_opt_3_loops.dat",force=True,quiet=True,numLoops=numLoops)
#    ocp.interpolateInitialGuess("data/crosswind_opt_3_loops_mechanical.dat",force=True,quiet=True,numLoops=1)
#    ocp.interpolateInitialGuess("data/crosswind_opt.dat",force=True,quiet=True,numLoops=numLoops)

    traj = ocp.solve()
#    from rawe.collocation import trajectory
#    traj = trajectory.TrajectoryPlotter(ocp,numpy.array(ocp._guess.vectorize()))


    print "num loops: "+str(numLoops)
    print "optimizing "+powerType
    print "optimal mechanical power: "+str(traj.lookup('mechanical_energy',-1)/traj.lookup('endTime'))
    print "optimal electrical power: "+str(traj.lookup('electrical_energy',-1)/traj.lookup('endTime'))
    print "endTime: "+str(traj.lookup('endTime'))

    traj.save("data/crosswind_opt_"+powerType+".dat")
#    traj.save("data/crosswind_opt_3_loops.dat")
#    traj.save("data/crosswind_opt_3_loops_electrical.dat")

    def printBoundsFeedback():
        xOpt = traj.dvMap.vectorize()
        lbx = ocp.solver.input('lbx')
        ubx = ocp.solver.input('ubx')
        ocp._bounds.printBoundsFeedback(xOpt,lbx,ubx,reportThreshold=0)
    printBoundsFeedback()

    def printConstraintsViolations():
        lbg = ocp.solver.input('lbg')
        ubg = ocp.solver.input('ubg')
        ocp._gfcn.setInput(traj.getDvs(),0)
        ocp._gfcn.evaluate()
        g = ocp._gfcn.output()
        ocp._constraints.printViolations(g,lbg,ubg,reportThreshold=0)#1e-9)
    printConstraintsViolations()

    # Plot the results
    def plotResults():
#        traj.subplot(['x','y','z'])
#        traj.subplot(['dx','dy','dz'])
#        traj.subplot([['aileron','elevator'],['daileron','delevator']],title='control surfaces')
        traj.subplot([['dddr'],['daileron','delevator']],title='control surfaces')
        traj.subplot(['wind_at_altitude','dr'],title='')
#        traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed',title='airspeed')
        traj.subplot([['alpha_deg','alphaTail_deg'],['beta_deg','betaTail_deg']])
        traj.subplot([['cL'],['cD','cD_tether'],['L_over_D','L_over_D_with_tether']],title='')
#        traj.subplot([['winch_power'], ['tether_tension'],['accel_g','accel_without_gravity_g']])
        traj.subplot([['rpm'],['dr']])
        traj.subplot([['tether_tension'],['torque']])
        traj.plot(['mechanical_winch_power', 'electrical_winch_power'])
        traj.plot('r')
#        traj.subplot([['ddx','ddy','ddz'],['accel','accel without gravity']])
#        traj.plot(["loyds_limit","loyds_limit_exact","neg_winch_power"])
#        traj.plot(["loyd's limit","-(winch power)"],title='')
        traj.subplot([['daileronCost'],['delevatorCost'],['dddrCost']])
        traj.subplot(['r','dr','ddr','dddr'])
#        traj.subplot(['w_bn_b_x','w_bn_b_y','w_bn_b_z'])
#        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
#        traj.plot('line_angle_deg')
#        traj.plot('quadrature_energy')
#        traj.subplot(['energy','quadrature_energy'])
#        traj.plot(['energy','quadrature_energy'])
#        traj.plot('nu')

        plt.show()
    plotResults()


    def plotPaper():
        plt.figure()
        plt.subplot(221)
        traj._plot('cL','',showLegend=False)
        plt.legend(['$C_L$'])
        plt.ylim([0,1.7])

        plt.subplot(223)
        traj._plot('L_over_D','',showLegend=False)
        plt.legend(['$L/D$'])
        plt.ylim([0,5.5])

        plt.subplot(222)
        traj._plot('wind_at_altitude','',showLegend=False)
        plt.ylabel('[m/s]')
        plt.ylim([7,9.2])
        plt.legend(['wind at altitude'])
        plt.subplot(224)
        traj._plot('dr','',showLegend=False)
        plt.ylabel('[m/s]')
        plt.ylim([-10,9])
        plt.legend(['$\dot{l}$'])

#        traj.subplot(['cL','L/D','dr'],title='')
#        traj.plot(["loyd's limit","loyd's limit (exact)","-(winch power)"])
#        traj.plot(["loyd's limit","-(winch power)"],title='')
        plt.figure()
        traj._plot("loyds_limit",'',showLegend=False,style=':')
        traj._plot("neg_winch_power",'',showLegend=False)
        plt.legend(["Loyd's limit","winch power"])
        plt.ylabel('power [W]')
        plt.ylim([-400,1100])
        plt.grid()


        plt.show()
#    plotPaper()
