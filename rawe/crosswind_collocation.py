import casadi as C
import matplotlib.pyplot as plt
import numpy
from numpy import pi
import zmq

from collocation import Coll,trajectory
from config import readConfig
import kiteutils
import kite_pb2
import kiteproto
import models

def setupOcp(dae,conf,publisher,nk=50,nicp=1,deg=4,collPoly='RADAU'):
    ocp = Coll(dae, nk=nk,nicp=nicp,deg=deg,collPoly=collPoly)
    
    print "setting up collocation..."
    ocp.setupCollocation(ocp.lookup('endTime'))

    print "moar setting up ocp..."
    
    # constrain invariants
    def constrainInvariantErrs():
        dcm = ocp.lookup('dcm',timestep=0)
        err = C.mul(dcm.T,dcm)
        ocp.constrain( C.veccat([err[0,0] - 1, err[1,1]-1, err[2,2] - 1, err[0,1], err[0,2], err[1,2]]), '==', 0, tag=
                       ("initial dcm orthonormal",None))
        ocp.constrain(ocp.lookup('c',timestep=0), '==', 0, tag=('c(0)==0',None))
        ocp.constrain(ocp.lookup('cdot',timestep=0), '==', 0, tag=('cdot(0)==0',None))
    constrainInvariantErrs()

    # constrain line angle
    for k in range(0,nk+1):
        ocp.constrain(kiteutils.getCosLineAngle(ocp,k),'>=',C.cos(55*pi/180), tag=('line angle',k))

    # constrain airspeed
    def constrainAirspeedAlphaBeta():
        for k in range(0,nk):
            ocp.constrain(ocp.lookup('airspeed',timestep=k), '>=', 10, tag=('airspeed',k))
            ocp.constrainBnds(ocp.lookup('alpha(deg)',timestep=k), (-10,30), tag=('alpha(deg)',k))
            ocp.constrainBnds(ocp.lookup('beta(deg)', timestep=k), (-10,10), tag=('beta(deg)',k))
    constrainAirspeedAlphaBeta()
    def constrainCl():
        for k in range(0,nk):
            ocp.constrain(ocp.lookup('cL',timestep=k), '<=', 2.0, tag=('cL',k))
    constrainCl()

    # constrain tether force
    for k in range(nk):
        ocp.constrain( ocp.lookup('tether tension',timestep=k,degIdx=1), '>=', 0, tag=('tether tension',k))
        ocp.constrain( ocp.lookup('tether tension',timestep=k,degIdx=ocp.deg), '>=', 0, tag=('tether tension',k))

    # make it periodic
    for name in [ "y","z",
                  "dy","dz",
                  "w1","w2","w3",
                  "r","dr",
                  'aileron','elevator'
                  ]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1), tag=('periodic diff state \"'+name+'"',None))

    # periodic attitude
    kiteutils.periodicDcm(ocp)

    # bounds
    ocp.bound('aileron',(-0.04,0.04))
    ocp.bound('elevator',(-0.1,0.1))
    ocp.bound('daileron',(-2.0,2.0))
    ocp.bound('delevator',(-2.0,2.0))

    ocp.bound('x',(-200,200))
    ocp.bound('y',(-200,200))
    if 'minAltitude' in conf:
        ocp.bound('z',(conf['minAltitude'],200))
    else:
        ocp.bound('z',(0.5,200))
    ocp.bound('r',(1,200))
    ocp.bound('dr',(-30,30))
    ocp.bound('ddr',(-500,500))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-70,70))

    for w in ['w1','w2','w3']:
        ocp.bound(w,(-4*pi,4*pi))

    ocp.bound('endTime',(0.5,20))
    ocp.bound('w0',(10,10))

    # boundary conditions
#    ocp.bound('y',(0,0),timestep=0,quiet=True)

    # guesses
    ocp.guess('endTime',5.4)
    ocp.guess('w0',10)

    # objective function
    obj = 0
    for k in range(nk):
        # control regularization
        ddr = ocp.lookup('ddr',timestep=k)
        daileron = ocp.lookup('daileron',timestep=k)
        delevator = ocp.lookup('delevator',timestep=k)
        
        daileronSigma = 0.1
        delevatorSigma = 0.1
        ddrSigma = 1.0
        
        ailObj = daileron*daileron / (daileronSigma*daileronSigma)
        eleObj = delevator*delevator / (delevatorSigma*delevatorSigma)
        winchObj = ddr*ddr / (ddrSigma*ddrSigma)
        
        obj += ailObj + eleObj + winchObj
    ocp.setObjective( 1e0*obj/nk + ocp.lookup('energy',timestep=-1)/ocp.lookup('endTime') )
    ocp.setQuadratureDdt('quadrature energy', 'winch power')

    return ocp


if __name__=='__main__':
    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    print "reading config..."
    conf = readConfig('config.ini','configspec.ini')
    
    print "creating model..."
    dae = models.crosswind(conf,extraParams=['endTime'])

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,publisher,nk=50)

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
                                                                 lineAlpha=0.2,
                                                                 zeroDelta=True) )
            mc = kite_pb2.MultiCarousel()
            mc.css.extend(list(kiteProtos))

            mc.messages.append("w0: "+str(traj.lookup('w0')))
            mc.messages.append("iter: "+str(self.iter))
            mc.messages.append("endTime: "+str(traj.lookup('endTime')))
            mc.messages.append("average power: "+str(traj.lookup('quadrature energy',timestep=-1)/traj.lookup('endTime'))+" W")

            # bounds feedback
#            lbx = ocp.solver.input(C.NLP_LBX)
#            ubx = ocp.solver.input(C.NLP_UBX)
#            ocp._bounds.printBoundsFeedback(xOpt,lbx,ubx,reportThreshold=0)
            
            publisher.send_multipart(["multi-carousel", mc.SerializeToString()])


    # solver
    solverOptions = [ ("expand_f",True)
                    , ("expand_g",True)
                    , ("generate_hessian",True)
#                     ,("qp_solver",C.NLPQPSolver)
#                     ,("qp_solver_options",{'nlp_solver': C.IpoptSolver, "nlp_solver_options":{"linear_solver":"ma57"}})
                    , ("linear_solver","ma57")
                    , ("max_iter",1000)
                    , ("tol",1e-8)
#                    , ("Timeout", 1e6)
#                    , ("UserHM", True)
#                    , ("ScaleConIter",True)
#                    , ("ScaledFD",True)
#                    , ("ScaledKKT",True)
#                    , ("ScaledObj",True)
#                    , ("ScaledQP",True)
                    ]
    
    print "setting up solver..."
    ocp.setupSolver( solverOpts=solverOptions,
                     callback=MyCallback() )

    ocp.interpolateInitialGuess("data/crosswind_homotopy.dat",force=True,quiet=True)
    traj = ocp.solve()

    
    print "optimal power: "+str(traj.lookup('quadrature energy',-1)/traj.lookup('endTime'))
    print "endTime: "+str(traj.lookup('endTime'))

    traj.save("data/crosswind_opt.dat")

    def printBoundsFeedback():
        xOpt = traj.dvMap.vectorize()
        lbx = ocp.solver.input(C.NLP_LBX)
        ubx = ocp.solver.input(C.NLP_UBX)
        ocp._bounds.printBoundsFeedback(xOpt,lbx,ubx,reportThreshold=0)
    printBoundsFeedback()

    lbg = ocp.solver.input(C.NLP_LBG)
    ubg = ocp.solver.input(C.NLP_UBG)
    ocp._gfcn.setInput(traj.getDvs(),0)
    ocp._gfcn.evaluate()
    g = ocp._gfcn.output()
    
    ocp._constraints.printViolations(g,lbg,ubg,reportThreshold=0)
    # Plot the results
    def plotResults():
        traj.subplot(['x','y','z'])
        traj.subplot(['dx','dy','dz'])
        traj.subplot([['aileron','elevator'],['daileron','delevator']],title='control surfaces')
        traj.subplot(['r','dr','ddr'])
        traj.subplot(['wind at altitude','dr','dx'])
        traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed',title='airspeed')
        traj.subplot([['alpha(deg)','alphaTail(deg)'],['beta(deg)','betaTail(deg)']])
        traj.subplot(['cL','cD','L/D'])
        traj.subplot(['winch power', 'tether tension'])
        traj.plot('energy')
        traj.subplot(['w1','w2','w3'])
        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        traj.plot('quadrature energy')
#        traj.subplot(['energy','quadrature energy'])
#        traj.plot(['energy','quadrature energy'])
        
        plt.show()
    plotResults()
