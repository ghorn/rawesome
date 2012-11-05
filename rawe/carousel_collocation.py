import copy
import casadi as C
import matplotlib.pyplot as plt
import numpy
from numpy import pi
import zmq
import pickle

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

    # constraint line angle
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
            ocp.constrain(airspeed,'>=',5)
            ocp.constrainBnds(alphaDeg,(-5,10))
            ocp.constrainBnds(betaDeg,(-10,10))
    constrainAirspeedAlphaBeta()

    # make it periodic
    for name in [ "y","z",
                  "dy","dz",
                  "w1","w2","w3",
                  "e11","e22","e33",
                  "ddelta",
                  "aileron","elevator",
                  "r","dr"]:
        ocp.constrain(ocp.lookup(name,timestep=0),'==',ocp.lookup(name,timestep=-1))

    # periodic attitude
#    kiteutils.periodicEulers(ocp)
#    kiteutils.periodicOrthonormalizedDcm(ocp)

    # bounds
    ocp.bound('aileron',(-0.04,0.04))
    ocp.bound('elevator',(-0.1,0.1))

    ocp.bound('daileron',(-2,2))
    ocp.bound('delevator',(-2,2))

    ocp.bound('x',(0.1,1000))
    ocp.bound('y',(-100,100))
    ocp.bound('z',(-0.5,7))
    ocp.bound('r',(0.5,10))
    ocp.bound('dr',(-10,10))
    ocp.bound('ddr',(0,0))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-70,70))

    for w in ['w1','w2','w3']:
        ocp.bound(w,(-4*pi,4*pi))

    ocp.bound('delta',(-0.01,1.01*2*pi))
    ocp.bound('ddelta',(-pi/8,8*pi))
    ocp.bound('tc',(-1000,1000))
    ocp.bound('endTime',(0.5,4.0))
    ocp.bound('w0',(10,10))
    ocp.bound('energy',(-1e6,1e6))

    # boundary conditions
    ocp.bound('delta',(0,0),timestep=0)
    ocp.bound('delta',(2*pi,2*pi),timestep=-1)
    ocp.bound('energy',(0,0),timestep=0,quiet=True)
    
    # objective function
    obj = 0
    for k in range(nk):
        u = ocp.uVec(k)
        ddr = ocp.lookup('ddr',timestep=k)
        tc = ocp.lookup('tc',timestep=k)
        aileron = ocp.lookup('aileron',timestep=k)
        elevator = ocp.lookup('elevator',timestep=k)
        
        aileronSigma = 0.1
        elevatorSigma = 0.1
        torqueSigma = 1000.0
        ddrSigma = 5.0
        
#        tc = tc - 390

        ailObj = aileron*aileron / (aileronSigma*aileronSigma)
        eleObj = elevator*elevator / (elevatorSigma*elevatorSigma)
        winchObj = ddr*ddr / (ddrSigma*ddrSigma)
        torqueObj = tc*tc / (torqueSigma*torqueSigma)
        
        obj += ailObj + eleObj + winchObj + torqueObj
    ocp.setObjective( C.sumAll(obj) )

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
            for k in range(0,nk):
                j = nicp*(deg+1)*k
                kiteProtos.append( kiteproto.toKiteProto(C.DMatrix(opt['x'][:,j]),
                                                         C.DMatrix(opt['u'][:,j]),
                                                         C.DMatrix(opt['p']),
                                                         conf['kite']['zt'],
                                                         conf['carousel']['rArm']) )
#            kiteProtos = [kiteproto.toKiteProto(C.DMatrix(opt['x'][:,k]),C.DMatrix(opt['u'][:,k]),C.DMatrix(opt['p']), conf['kite']['zt'], conf['carousel']['rArm']) for k in range(opt['x'].shape[1])]
            
            mc = kite_pb2.MultiCarousel()
            mc.css.extend(list(kiteProtos+oldKites))
            
            mc.messages.append("endTime: "+str(xup['endTime']))
            mc.messages.append("w0: "+str(xup['w0']))
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

    ocp.guess('tc',0)
    ocp.guess('endTime',1.5)
    ocp.guess('daileron',0)
    ocp.guess('delevator',0)

    ocp.guess('ddr',0)
    ocp.guess('w0',10)

    ocp.setupCollocation(ocp.lookup('endTime'))
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
    ocp = setupOcp(dae,conf,publisher,nk=50)

    xOpt = None
    for w0 in [10]:
        ocp.bound('w0',(w0,w0),force=True)
        opt = ocp.solve(xInit=xOpt)
        xup = opt['vardict']
        xOpt = opt['X_OPT']
        
        for k in range(0,ocp.nk):
            j = ocp.nicp*(ocp.deg+1)*k
            oldKites.append( kiteproto.toKiteProto(C.DMatrix(opt['x'][:,j]),C.DMatrix(opt['u'][:,j]),C.DMatrix(opt['p']), conf['kite']['zt'], conf['carousel']['rArm']) )

    print "saving optimal trajectory"
    f=open("data/carousel_opt.dat",'w')
    pickle.dump(opt,f)
    f.close()

    print "optimal power: "+str(opt['vardict']['energy'][-1]/opt['vardict']['endTime'])
    # Plot the results
    ocp.plot(['x','y','z'],opt)
    ocp.plot(['aileron','elevator'],opt,title='control surface inputs')
    ocp.plot(['tc'],opt,title='motor inputs (tc)')
    ocp.plot(['ddr'],opt,title='winch accel (ddr)')
    ocp.subplot(['c','cdot','cddot'],opt,title="invariants")
    ocp.plot('airspeed',opt)
    ocp.subplot([['alpha(deg)','alphaTail(deg)'],['beta(deg)','betaTail(deg)']],opt)
    ocp.subplot(['cL','cD','L/D'],opt)
    ocp.plot('motor power',opt)
#    ocp.plot('winch power',opt)
    ocp.plot('energy',opt)
    ocp.subplot(['e11','e12','e13',
                 'e21','e22','e23',
                 'e31','e32','e33'],opt)
    plt.show()
