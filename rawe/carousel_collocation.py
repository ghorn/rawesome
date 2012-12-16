import copy
import casadi as C
import matplotlib.pyplot as plt
import numpy
from numpy import pi
import zmq
import pickle

import trajectory
from collocation import Coll,boundsFeedback
from config import readConfig
import kiteutils
import kite_pb2
import kiteproto
import models

x0 = C.DMatrix( [ 1.154244772411
                , -0.103540608242
                , -0.347959211327
                , 0.124930983341
                , 0.991534857363
                , 0.035367725910
                , 0.316039689643
                , -0.073559821379
                , 0.945889986864
                , 0.940484536806
                , -0.106993361072
                , -0.322554269411
                , 0.000000000000
                , 0.000000000000
                , 0.000000000000
                , 0.137035790811
                , 3.664945343102
                , -1.249768772258
                , 0.000000000000
                , 3.874600000000
                , 0.0
                ])
x0=C.veccat([x0,C.sqrt(C.sumAll(x0[0:2]*x0[0:2])),0,0,0])

oldKites = []

def setupOcp(dae,conf,publisher,nk=50,nicp=1,deg=4):
    ocp = Coll(dae, nk=nk,nicp=nicp,deg=deg)
    ocp.setupCollocation(ocp.lookup('endTime'))
                   
    # constrain invariants
    def constrainInvariantErrs():
        dcm = ocp.lookup('dcm',timestep=0)
        err = C.mul(dcm.T,dcm)
        ocp.constrain( C.veccat([err[0,0] - 1, err[1,1]-1, err[2,2] - 1, err[0,1], err[0,2], err[1,2]]), '==', 0)
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0)
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0)
    constrainInvariantErrs()

    # constraint line angle
    for k in range(0,nk+1):
        ocp.constrain(kiteutils.getCosLineAngle(ocp,k),'>=',C.cos(55*pi/180))
        
    # constrain airspeed
    def constrainAirspeedAlphaBeta():
        for k in range(0,nk):
            ocp.constrain(ocp.lookup('airspeed',timestep=k), '>=', 5)
            ocp.constrainBnds(ocp.lookup('alpha(deg)',timestep=k), (-5,10))
            ocp.constrainBnds(ocp.lookup('beta(deg)', timestep=k), (-10,10))
    constrainAirspeedAlphaBeta()

    # constrain tether force
    for k in range(nk):
        ocp.constrain( ocp.lookup('tether tension',timestep=k,degIdx=1), '>=', 0)
        ocp.constrain( ocp.lookup('tether tension',timestep=k,degIdx=ocp.deg), '>=', 0)

    # make it periodic
    for name in [ "y","z",
                  "dy","dz",
                  "w1","w2","w3",
                  "ddelta",
                  "aileron","elevator",
                  "r","dr"]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1))

    # periodic attitude
#    kiteutils.periodicEulers(ocp)
#    kiteutils.periodicOrthonormalizedDcm(ocp)
    kiteutils.periodicDcm(ocp)

    # bounds
    ocp.bound('aileron',(-0.04,0.04))
    ocp.bound('elevator',(-0.1,0.1))
    ocp.bound('daileron',(-2,2))
    ocp.bound('delevator',(-2,2))

    ocp.bound('x',(0.1,1000))
    ocp.bound('y',(-100,100))
    ocp.bound('z',(-0.5,7))
    ocp.bound('r',(4,4))
    ocp.bound('dr',(-10,10))
    ocp.bound('ddr',(-2.5,2.5))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-70,70))

    for w in ['w1','w2','w3']:
        ocp.bound(w,(-4*pi,4*pi))

    ocp.bound('delta',(-0.01,1.01*2*pi))
    ocp.bound('ddelta',(-pi/8,8*pi))
    ocp.bound('motor torque',(-1000,1000))
    ocp.bound('endTime',(0.5,7.0))
    ocp.bound('w0',(10,10))
    ocp.bound('energy',(-1e6,1e6))

    # boundary conditions
    ocp.bound('delta',(0,0),timestep=0)
    ocp.bound('delta',(2*pi,2*pi),timestep=-1)
    ocp.bound('energy',(0,0),timestep=0,quiet=True)
    
    # objective function
    obj = 0
    for k in range(nk):
        # control regularization
        ddr = ocp.lookup('ddr',timestep=k)
        tc = ocp.lookup('motor torque',timestep=k)
        daileron = ocp.lookup('daileron',timestep=k)
        delevator = ocp.lookup('delevator',timestep=k)
        
        daileronSigma = 0.1
        delevatorSigma = 0.1
        ddrSigma = 5.0
        torqueSigma = 1.0
        
#        tc = tc - 390

        ailObj = daileron*daileron / (daileronSigma*daileronSigma)
        eleObj = delevator*delevator / (delevatorSigma*delevatorSigma)
        winchObj = ddr*ddr / (ddrSigma*ddrSigma)
        torqueObj = tc*tc / (torqueSigma*torqueSigma)
        
        obj += ailObj + eleObj + winchObj + torqueObj
    ocp.setObjective( obj/nk )

    # callback function
    class MyCallback:
        def __init__(self):
            self.iter = 0 
        def __call__(self,f,*args):
            self.iter = self.iter + 1
            xOpt = numpy.array(f.input(C.NLP_X_OPT))

            traj = trajectory.Trajectory(ocp,xOpt)
            
            kiteProtos = []
            for k in range(0,ocp.nk):
                for nicpIdx in range(0,ocp.nicp):
                    for j in [0]:
#                    for j in range(ocp.deg+1):
                        kiteProtos.append( kiteproto.toKiteProto(C.DMatrix(traj.dvMap.xVec(k,nicpIdx=nicpIdx,degIdx=j)),
                                                                 C.DMatrix(traj.dvMap.uVec(k)),
                                                                 C.DMatrix(traj.dvMap.pVec()),
                                                                 conf['kite']['zt'],
                                                                 conf['carousel']['rArm'],
                                                                 lineAlpha=0.2) )
            mc = kite_pb2.MultiCarousel()
            mc.css.extend(list(kiteProtos+oldKites))
            
            mc.messages.append("w0: "+str(traj.lookup('w0')))
            mc.messages.append("endTime: "+str(traj.lookup('endTime')))
            mc.messages.append("iter: "+str(self.iter))

#            # bounds feedback
#            lbx = ocp.solver.input(C.NLP_LBX)
#            ubx = ocp.solver.input(C.NLP_UBX)
#
#            violations = boundsFeedback(xOpt,lbx,ubx,ocp.bndtags)
#            for name in violations:
#                print "violation!: "+name+": "+str(violations[name])
#                mc.messages.append(name+": "+str(violations[name]))
            
            publisher.send_multipart(["multi-carousel", mc.SerializeToString()])


    # solver
    solverOptions = [ ("expand_f",True)
                    , ("expand_g",True)
                    , ("generate_hessian",True)
#                     ,("qp_solver",C.NLPQPSolver)
#                     ,("qp_solver_options",{'nlp_solver': C.IpoptSolver, "nlp_solver_options":{"linear_solver":"ma57"}})
                    , ("linear_solver","ma57")
                    , ("max_iter",1000)
                    , ("tol",1e-9)
#                    , ("Timeout", 1e6)
#                    , ("UserHM", True)
#                    , ("ScaleConIter",True)
#                    , ("ScaledFD",True)
#                    , ("ScaledKKT",True)
#                    , ("ScaledObj",True)
#                    , ("ScaledQP",True)
                    ]
    
    # initial guess
    ocp.guessX(x0)
    for k in range(0,nk+1):
        val = 2.0*pi*k/nk
        ocp.guess('delta',val,timestep=k,quiet=True)

    ocp.guess('motor torque',0)
    ocp.guess('endTime',1.5)
    ocp.guess('daileron',0)
    ocp.guess('delevator',0)

    ocp.guess('ddr',0)
    ocp.guess('w0',10)

    ocp.setupSolver( solverOpts=solverOptions,
                     callback=MyCallback() )

    return ocp


if __name__=='__main__':
    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    print "reading config..."
    conf = readConfig('config.ini','configspec.ini')
    
    print "creating model..."
    dae = models.carousel(conf,extraParams=['endTime'])

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,publisher,nk=30)

    ocp.interpolateInitialGuess("data/carousel_opt.dat",force=True,quiet=True)

    for w0 in [10]:
        ocp.bound('w0',(w0,w0),force=True)
        traj = ocp.solve()
        
        for k in range(0,ocp.nk):
            for nicpIdx in range(0,ocp.nicp):
                for j in [0]:
#                for j in range(ocp.deg+1):
                    oldKites.append( kiteproto.toKiteProto(C.DMatrix(traj.dvMap.xVec(k,nicpIdx=nicpIdx,degIdx=j)),
                                                           C.DMatrix(traj.dvMap.uVec(k)),
                                                           C.DMatrix(traj.dvMap.pVec()),
                                                           conf['kite']['zt'],
                                                           conf['carousel']['rArm'],
                                                           lineAlpha=0.2) )

    print "saving optimal trajectory"
    traj.save("data/carousel_opt.dat")

    print "optimal power: "+str(traj.lookup('energy',timestep=-1)/traj.lookup('endTime'))
    
    # Plot the results
    def plotResults():
        traj.plot(['x','y','z'])
        traj.plot(['aileron','elevator'],title='control surface inputs')
        traj.plot(['ddr'],title='winch accel (ddr)')
        traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed')
        traj.subplot([['alpha(deg)','alphaTail(deg)'],['beta(deg)','betaTail(deg)']])
        traj.subplot(['cL','cD','L/D'])
        traj.subplot(['motor torque','motor power'])
        traj.subplot(['winch power','tether tension'])
        traj.plot('energy')
        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        plt.show()
    plotResults()
