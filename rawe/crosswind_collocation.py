import casadi as C
import matplotlib.pyplot as plt
import numpy
from numpy import pi
import zmq
import pickle
from trajectory import Trajectory

from collocation import Coll,boundsFeedback
from config import readConfig
import kiteutils
import kite_pb2
import kiteproto
import models

def setupOcp(dae,conf,publisher,nk=50,nicp=1,deg=4):
    ocp = Coll(dae, nk=nk,nicp=nicp,deg=deg)
    
    # constrain invariants
    def invariantErrs():
        dcm = C.horzcat( [ C.veccat([dae.x('e11'), dae.x('e21'), dae.x('e31')])
                         , C.veccat([dae.x('e12'), dae.x('e22'), dae.x('e32')])
                         , C.veccat([dae.x('e13'), dae.x('e23'), dae.x('e33')])
                         ] ).trans()
        err = C.mul(dcm.trans(), dcm)
        dcmErr = C.veccat([ err[0,0]-1, err[1,1]-1, err[2,2]-1, err[0,1], err[0,2], err[1,2] ])
        f = C.SXFunction( [dae.xVec(),dae.uVec(),dae.pVec()]
                        , [dae.output('c'),dae.output('cdot'),dcmErr]
                        )
        f.setOption('name','invariant errors')
        f.init()
        return f

    [c0,cdot0,dcmError0] = invariantErrs().call([ocp.xVec(0),ocp.uVec(0),ocp.pVec()])
    ocp.constrain(c0,'==',0)
    ocp.constrain(cdot0,'==',0)
    ocp.constrain(dcmError0,'==',0)

    # constrain line angle
    for k in range(0,nk+1):
        ocp.constrain(kiteutils.getCosLineAngle(ocp,k),'>=',C.cos(55*pi/180))

    # constrain airspeed
    def constrainAirspeedAlphaBeta():
        f = C.SXFunction( [dae.xVec(),dae.uVec(),dae.pVec()]
                        , [dae.output('airspeed'),dae.output('alpha(deg)'),dae.output('beta(deg)')]
                        )
        f.setOption('name','airspeed/alpha/beta')
        f.init()

        for k in range(0,nk):
            [airspeed,alphaDeg,betaDeg] = f.call([ocp.xVec(k),ocp.uVec(k),ocp.pVec()])
            ocp.constrain(airspeed,'>=',10)
            ocp.constrainBnds(alphaDeg,(-5,10))
            ocp.constrainBnds(betaDeg,(-10,10))
    constrainAirspeedAlphaBeta()

    # constrain tether force
    for k in range(nk):
        degIdx = 1# in range(1,deg+1):
        r  = ocp.lookup('r', timestep=k,degIdx=degIdx)
        nu = ocp.lookup('nu',timestep=k,degIdx=degIdx)
        ocp.constrain(0,'<=',r*nu)

        degIdx = deg# in range(1,deg+1):
        r  = ocp.lookup('r', timestep=k,degIdx=degIdx)
        nu = ocp.lookup('nu',timestep=k,degIdx=degIdx)
        ocp.constrain(0,'<=',r*nu)

    # make it periodic
    for name in [ "y","z",
                  "dy","dz",
                  "w1","w2","w3",
                  "e11","e22","e33",
                  "r","dr",
                  'aileron','elevator'
                  ]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1))

    # make sure it doesn't find transposed-periodic DCM
    # sign(eij(beginning) == sign(eij(end)) <--> eij(beginning)*eij(end) >= 0
    for name in ['e12','e13','e23']:
        ocp.constrain(ocp.lookup(name,timestep=0)*ocp.lookup(name,timestep=-1),'>=',0)

    # periodic attitude
#    kiteutils.periodicEulers(ocp)
#    kiteutils.periodicOrthonormalizedDcm(ocp)

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
    ocp.bound('r',(1,100))
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
    ocp.bound('energy',(-1e6,1e6))

    # boundary conditions
    ocp.bound('energy',(0,0),timestep=0,quiet=True)
#    ocp.bound('y',(0,0),timestep=0,quiet=True)
    
    
    # initial guess
#    ocp.guessX(x0)
#    for k in range(0,nk+1):
#        val = 2.0*pi*k/nk
#        ocp.guess('delta',val,timestep=k,quiet=True)
#
#    ocp.guess('aileron',0)
#    ocp.guess('elevator',0)
#    ocp.guess('tc',0)
    ocp.guess('endTime',5.4)
#
#    ocp.guess('ddr',0)
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

    return ocp


if __name__=='__main__':
    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    print "reading config..."
    conf = readConfig('config.ini','configspec.ini')
    
    print "creating model..."
    dae = models.crosswind(conf,extraParams=['endTime'])

    # load initial guess
    print "loading initial guess data..."
    f = open('data/crosswind_guess.txt','r')
    xutraj = []
    for line in f:
        xutraj.append([float(x) for x in line.strip('[]\n').split(',')])
    f.close()
    xutraj = numpy.array(xutraj)
    # remove delta/ddelta/tc
    xutraj = numpy.hstack((xutraj[:,:23], xutraj[:,24:]))
    xutraj = numpy.hstack((xutraj[:,:18], xutraj[:,20:]))

    # add daileron/delevator
    xutraj = numpy.hstack((xutraj, 0*xutraj[:,:2]))

    # make it go around n times
    nTurns = 1
    xutraj = numpy.vstack(tuple([xutraj]*nTurns))
    
    print "setting up ocp..."
    ocp = setupOcp(dae,conf,publisher,nk=100)

    # callback function
    class MyCallback:
        def __init__(self):
            self.iter = 0 
        def __call__(self,f,*args):
            self.iter = self.iter + 1
            xOpt = numpy.array(f.input(C.NLP_X_OPT))

            opt = ocp.devectorize(xOpt)
            xup = opt['vardict']
            
            kiteProtos = []
            for k in range(0,ocp.nk):
                j = ocp.nicp*(ocp.deg+1)*k
                kiteProtos.append( kiteproto.toKiteProto(C.DMatrix(opt['x'][:,j]),
                                                         C.DMatrix(opt['u'][:,j]),
                                                         C.DMatrix(opt['p']),
                                                         conf['kite']['zt'],
                                                         conf['carousel']['rArm'],
                                                         lineAlpha=0.2,
                                                         zeroDelta=True) )
#            kiteProtos = [kiteproto.toKiteProto(C.DMatrix(opt['x'][:,k]),C.DMatrix(opt['u'][:,k]),C.DMatrix(opt['p']), conf['kite']['zt'], conf['carousel']['rArm'], zeroDelta=True) for k in range(opt['x'].shape[1])]
            
            mc = kite_pb2.MultiCarousel()
            mc.css.extend(list(kiteProtos))
            
            mc.messages.append("w0: "+str(xup['w0']))
            mc.messages.append("iter: "+str(self.iter))
            mc.messages.append("endTime: "+str(xup['endTime']))
            mc.messages.append("average power: "+str(xup['energy'][-1]/xup['endTime'])+" W")

            # bounds feedback
#            lbx = ocp.solver.input(C.NLP_LBX)
#            ubx = ocp.solver.input(C.NLP_UBX)
#            violations = boundsFeedback(xOpt,lbx,ubx,ocp.bndtags,tolerance=1e-9)
#            for name in violations:
#                violmsg = "violation!: "+name+": "+str(violations[name])
#                mc.messages.append(violmsg)
            
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
    
    print "setting up collocation..."
    ocp.setupCollocation(ocp.lookup('endTime'))
    print "setting up solver..."
    ocp.setupSolver( solverOpts=solverOptions,
                     callback=MyCallback() )

    print "interpolating initial guess..."
    xuguess = numpy.array([numpy.interp(numpy.linspace(0,1,ocp.nk+1), numpy.linspace(0,1,xutraj.shape[0]), xutraj[:,k]) for k in range(xutraj.shape[1])])
    for k in range(ocp.nk+1):
        ocp.guessX(xuguess[:len(ocp.dae.xNames()),k],timestep=k,quiet=True)
        if k < ocp.nk:
            ocp.guessU(xuguess[len(ocp.dae.xNames()):,k],timestep=k,quiet=True)
    ocp.guess('energy',0,quiet=True)
    
    print "loading optimal trajectory"
    f=open("data/crosswind_opt.dat",'r')
    opt = pickle.load(f)
    f.close()

#    opt = ocp.solve(xInit=opt['X_OPT'])
    opt = ocp.solve()
    traj = Trajectory(ocp,dvs=opt['X_OPT'])
    
    print "optimal power: "+str(opt['vardict']['energy'][-1]/opt['vardict']['endTime'])
    
    print "saving optimal trajectory"
    traj.save("data/crosswind_opt.dat")

    def printBoundsFeedback():
#        (_,lbx,ubx) = ocp._vectorizeBoundsAndGuess( ocp._parseBoundsAndGuess(False,False) )
        lbx = ocp.solver.input(C.NLP_LBX)
        ubx = ocp.solver.input(C.NLP_UBX)
        violations = boundsFeedback(opt['X_OPT'],lbx,ubx,ocp.bndtags,tolerance=0.5)
        for name in violations:
            print "violation!: "+name+": "+str(violations[name])
#    printBoundsFeedback()

    # Plot the results
    def plotResults():
        traj.subplot(['x','y','z'])
        traj.subplot(['dx','dy','dz'])
        traj.subplot([['aileron','elevator'],['daileron','delevator']],title='control surfaces')
        traj.subplot(['r','dr','ddr'])
        traj.subplot(['wind at altitude','dr','dx'])
        traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed')
        traj.subplot([['alpha(deg)','alphaTail(deg)'],['beta(deg)','betaTail(deg)']])
        traj.subplot(['cL','cD','L/D'])
        traj.subplot(['winch power', 'tether tension'])
        traj.plot('energy')
        traj.subplot(['w1','w2','w3'])
        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        plt.show()
    plotResults()
