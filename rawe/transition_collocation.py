import copy
import casadi as C
import matplotlib.pyplot as plt
import numpy
from numpy import pi
import zmq
from fourier_fit import FourierFit
import pickle
from trajectory import Trajectory

from collocation import Coll,boundsFeedback
from config import readConfig
import kiteutils
import kite_pb2
import kiteproto
import models

def toFourierKiteProto(fits,k,kiteAlpha=1,lineAlpha=1,dz=0,zt=0,rArm=0):
    cs = kite_pb2.CarouselState()

    cs.kiteXyz.x = fits['x'].Mx[k]
    cs.kiteXyz.y = fits['y'].Mx[k]
    cs.kiteXyz.z = fits['z'].Mx[k]+dz

    cs.kiteDcm.r11 = fits['e11'].Mx[k]
    cs.kiteDcm.r12 = fits['e12'].Mx[k]
    cs.kiteDcm.r13 = fits['e13'].Mx[k]

    cs.kiteDcm.r21 = fits['e21'].Mx[k]
    cs.kiteDcm.r22 = fits['e22'].Mx[k]
    cs.kiteDcm.r23 = fits['e23'].Mx[k]

    cs.kiteDcm.r31 = fits['e31'].Mx[k]
    cs.kiteDcm.r32 = fits['e32'].Mx[k]
    cs.kiteDcm.r33 = fits['e33'].Mx[k]

    if 'delta' in fits:
        cs.delta = fits['delta'].Mx[k]
    else:
        cs.delta = 0

    cs.rArm = rArm
    cs.zt = zt

    cs.kiteTransparency = kiteAlpha
    cs.lineTransparency = lineAlpha

    return cs


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
            ocp.constrain(airspeed,'>=',5)
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

    # bounds
    ocp.bound('aileron',(-0.04,0.04))
    ocp.bound('elevator',(-0.1,0.1))
    ocp.bound('daileron',(-2.0,2.0))
    ocp.bound('delevator',(-2.0,2.0))
    ocp.bound('motor torque',(-10000,10000))

    ocp.bound('x',(-200,200))
    ocp.bound('y',(-200,200))
    ocp.bound('z',(-0.5,200))
    ocp.bound('r',(1,31))
    ocp.bound('dr',(-30,30))
    ocp.bound('ddr',(-500,500))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        ocp.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        ocp.bound(d,(-70,70))

    for w in ['w1','w2','w3']:
        ocp.bound(w,(-4*pi,4*pi))

    ocp.bound('delta',(-12*pi,12*pi))
    ocp.bound('ddelta',(-12*pi,12*pi))
    ocp.bound('endTime',(1.5,15))
    ocp.bound('w0',(10,10))

    ocp.bound('phase0',(-12*pi,12*pi))
    ocp.bound('phaseF',(-12*pi,12*pi))

    ocp.bound('energy',(-1e6,1e6))
    
    # boundary conditions
    ocp.bound('energy',(0,0),timestep=0,quiet=True)
    ocp.bound('delta',(0,0),timestep=-1,quiet=True)
    ocp.bound('ddelta',(0,0),timestep=-1,quiet=True)

    def getFourierFit(filename,phase):
        f=open(filename,'r')
        fits0 = pickle.load(f)
        f.close()

        fourierBases = []
        for po in fits0['x'].polyOrder:
            fourierBases.append(phase**po)
        for co in fits0['x'].cosOrder:
            fourierBases.append(C.cos(co*phase))
        for so in fits0['x'].sinOrder:
            fourierBases.append(C.sin(so*phase))
        fourierFit = {}
        for name in fits0:
            fitCoeffs = fits0[name].fitcoeffs
            fourierFit[name] = sum([fc*fb for fc,fb in zip(fitCoeffs,fourierBases)])

        return (fourierFit,fits0)
        
    (startup,startupfits) = getFourierFit("data/carousel_opt_fourier.dat",ocp.lookup('phase0'))
    (crosswind,crosswindfits) = getFourierFit("data/crosswind_opt_fourier.dat",ocp.lookup('phaseF'))

#    def matchStartupEuler():
#        (yaw0,pitch0,roll0) = kiteutils.getEuler(ocp, 0)
#        ocp.constrain(yaw0,'==',startup['yaw'])
#        ocp.constrain(pitch0,'==',startup['pitch'])
#        ocp.constrain(roll0,'==',startup['roll'])
#    matchStartupEuler()
    
#    for name in ['e11','e22','e33','y','z','dy','dz','r','dr','w1','w2','w3','delta','ddelta']:
    for name in ['e11','e22','e33','y','z','dy','dz','r','dr','delta','ddelta']:
        ocp.constrain(ocp.lookup(name,timestep=0), '==', startup[name])
    
#    def matchFinalEuler():
#        (yawf,pitchf,rollf) = kiteutils.getEuler(ocp, -1)
#        ocp.constrain(yawf,'==',crosswind['yaw'])
#        ocp.constrain(pitchf,'==',crosswind['pitch'])
#        ocp.constrain(rollf,'==',crosswind['roll'])
#    matchFinalEuler()
#    for name in ['e11','e22','e33']:
#        ocp.constrain(ocp.lookup(name,timestep=-1), '==', crosswind[name])
    
#    for name in ['e11','e22','e33','y','z','dy','dz','r','dr','w1','w2','w3']:
    for name in ['e11','e22','e33','y','z','dy','dz','r','dr']:
        ocp.constrain(ocp.lookup(name,timestep=-1), '==', crosswind[name])

    # make sure it doesn't find transposed-periodic DCM
    # sign(eij(beginning) == sign(eij(end)) <--> eij(beginning)*eij(end) >= 0
    for name in ['e12','e13','e23']:
        ocp.constrain(ocp.lookup(name,timestep=0) *  startup[name],'>=',0)
        ocp.constrain(ocp.lookup(name,timestep=-1)*crosswind[name],'>=',0)

    # initial guess
    phase0Guess = -4.0*pi
    phaseFGuess = -0.4*2*pi
    
    namesF = ['x','y','z','dx','dy','dz','r','dr','w1','w2','w3','e11','e12','e13','e21','e22','e23','e31','e32','e33']
    names0 = namesF+['delta','ddelta']

    initfun = C.MXFunction([ocp.lookup('phase0')],[startup[name] for name in names0])
    initfun.init()
    for j in range(nk+1):
        alpha = float(j)/nk
        phase = phase0Guess*(1-alpha) + phaseFGuess*alpha
        initfun.setInput([phase])
        initfun.evaluate()
        initvals = dict([(name,float(initfun.output(k)[0,0])) for k,name in enumerate(names0)])
        for name in initvals:
            ocp.guess(name,float(initvals[name]),timestep=j)

#    finalfun = C.MXFunction([ocp.lookup('phaseF')],[crosswind[name] for name in namesF])
#    finalfun.init()
#    finalfun.setInput([phaseFGuess])
#    finalfun.evaluate()
#    finalvals = dict([(name,float(finalfun.output(k)[0,0])) for k,name in enumerate(namesF)])
#    initvals['delta'] = 0
#    initvals['ddelta'] = 0
#    finalvals['delta'] = 0
#    finalvals['ddelta'] = 0

#    for name in namesF:
#        for k in range(nk+1):
#            alpha = float(k)/nk
#            ocp.guess(name,(1-alpha)*initvals[name] + alpha*finalvals[name], timestep=k)
#    
#    for j,name in enumerate(namesF):
#        for k in range(nk+1):
#            alpha = float(k)/nk
#            finalfun.setInput([phaseFGuess])
#            finalval = float(finalfun.output(k)[0,0])) for k,name in enumerate(namesF)])
#            ocp.guess(name,(1-alpha)*initvals[name] + alpha*finalvals[name], timestep=k)
    
    ocp.guess('energy',0)
    ocp.guess('motor torque',0)
    ocp.guess('aileron',0)
    ocp.guess('elevator',0)
    ocp.guess('ddr',0)
    ocp.guess('phase0',phase0Guess)
    ocp.guess('phaseF',phaseFGuess)
    ocp.guess('endTime',8)
    ocp.guess('w0',10)
    
    # objective function
    obj = 0
    for k in range(nk):
        u = ocp.uVec(k)
        ddr = ocp.lookup('ddr',timestep=k)
        tc = ocp.lookup('motor torque',timestep=k)
        torqueSigma = 1000.0
        aileron = ocp.lookup('aileron',timestep=k)
        elevator = ocp.lookup('elevator',timestep=k)
        
        aileronSigma = 0.1
        elevatorSigma = 0.1
        torqueSigma = 1000.0
        ddrSigma = 5.0
        
        ailObj = aileron*aileron / (aileronSigma*aileronSigma)
        eleObj = elevator*elevator / (elevatorSigma*elevatorSigma)
        winchObj = ddr*ddr / (ddrSigma*ddrSigma)
        torqueObj = tc*tc / (torqueSigma*torqueSigma)
        
        obj += ailObj + eleObj + winchObj + torqueObj
        
    ocp.setObjective( 1e1*obj/nk + ocp.lookup('endTime') )

    oldKiteProtos = []
    for k in range(0,startupfits['x'].Mx.size):
        oldKiteProtos.append( toFourierKiteProto(startupfits,k,zt=conf['kite']['zt'], rArm=conf['carousel']['rArm'],kiteAlpha=0.3,lineAlpha=0,dz=0) )
    for k in range(0,crosswindfits['x'].Mx.size):
        oldKiteProtos.append( toFourierKiteProto(crosswindfits,k,zt=conf['kite']['zt'], rArm=conf['carousel']['rArm'],kiteAlpha=0.3,lineAlpha=0,dz=0) )

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
                                                         conf['carousel']['rArm'],
                                                         lineAlpha=0.2)
                                   )
#            kiteProtos = [kiteproto.toKiteProto(C.DMatrix(opt['x'][:,k]),C.DMatrix(opt['u'][:,k]),C.DMatrix(opt['p']), conf['kite']['zt'], conf['carousel']['rArm'], zeroDelta=True) for k in range(opt['x'].shape[1])]
            kiteProtos += oldKiteProtos
            mc = kite_pb2.MultiCarousel()
            mc.css.extend(list(kiteProtos))
            
            mc.messages.append("endTime: "+str(xup['endTime']))
            mc.messages.append("phase0: "+str(xup['phase0']/pi)+" * pi")
            mc.messages.append("phaseF: "+str(xup['phaseF']/pi)+" * pi")
            mc.messages.append("w0: "+str(xup['w0']))
            mc.messages.append("iter: "+str(self.iter))

#            # bounds feedback
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
                    , ("tol",1e-9)
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
    return ocp


if __name__=='__main__':
    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    print "reading config..."
    conf = readConfig('config.ini','configspec.ini')
    
    print "creating model..."
    dae = models.carousel(conf,extraParams=['endTime','phase0','phaseF'])

    # load initial guess
    print "loading initial guess data..."
    f = open('data/transition_guess.txt','r')
    xutraj = []
    for line in f:
        xutraj.append([float(x) for x in line.strip('[]\n').split(',')])
    f.close()
    xutraj = numpy.array(xutraj)
    xutraj = xutraj[220:-10,:]
    xutraj[:,18] = xutraj[:,18] - 6*pi
    # remove delta/ddelta/tc
#    xutraj = numpy.hstack((xutraj[:,:23], xutraj[:,24:]))
#    xutraj = numpy.hstack((xutraj[:,:18], xutraj[:,20:]))
    # add aileron/elevator
#    xutraj = numpy.hstack((xutraj[:,:23], xutraj[:,23:]))
    xutraj = numpy.hstack((xutraj[:,:23], 0*xutraj[:,:2], xutraj[:,23:]))

    print "setting up ocp..."
    ocp = setupOcp(dae,conf,publisher,nk=80)

    print "interpolating initial guess..."
    xuguess = numpy.array([numpy.interp(numpy.linspace(0,1,ocp.nk+1), numpy.linspace(0,1,xutraj.shape[0]), xutraj[:,k]) for k in range(xutraj.shape[1])])
    for k in range(ocp.nk+1):
        ocp.guessX(xuguess[:len(ocp.dae.xNames()),k],timestep=k,quiet=True)
        if k < ocp.nk:
            ocp.guessU(xuguess[len(ocp.dae.xNames()):,k],timestep=k,quiet=True)
    ocp.guess('energy',0,quiet=True)

    # load old initial guess
    f=open("data/transition_working_guess.dat",'r')
    workingGuess = pickle.load(f)
    f.close()
    opt = ocp.solve(xInit=workingGuess['X_OPT'])
    print "optimal power: "+str(opt['vardict']['energy'][-1]/opt['vardict']['endTime'])
    
    traj = Trajectory(ocp,dvs=opt['X_OPT'])
    print "saving optimal trajectory"
    traj.save("data/transition_opt.dat")
    
    # Plot the results
    def plotResults():
##        traj.plot(['x','y','z'])
        traj.subplot([['aileron','elevator'],['daileron','delevator']],title='control surfaces')
        traj.subplot(['r','dr','ddr'])
        traj.subplot(['c','cdot','cddot'],title="invariants")
        traj.plot('airspeed')
        traj.subplot([['alpha(deg)','alphaTail(deg)'],['beta(deg)','betaTail(deg)']])
        traj.subplot(['cL','cD','L/D'])
        traj.subplot(['winch power', 'tether tension'])
        traj.plot('motor torque')
        traj.plot('energy')
#        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        plt.show()
    plotResults()
