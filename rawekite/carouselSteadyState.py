import zmq
import casadi as C

from rawe.ocputils import Constraints

pi = C.pi

def getSteadyState(dae,conf,omega0,r0):
    # make steady state model
    g = Constraints()
    g.add(dae.getResidual(),'==',0,tag=('dae residual',None))
    def constrainInvariantErrs():
        dcm = dae['dcm']
        err = C.mul(dcm.T,dcm)
        g.add( C.veccat([err[0,0] - 1, err[1,1]-1, err[2,2] - 1, err[0,1], err[0,2], err[1,2]]), '==', 0, tag=
                       ("initial dcm orthonormal",None))
        g.add(dae['c'], '==', 0, tag=('c(0)==0',None))
        g.add(dae['cdot'], '==', 0, tag=('cdot(0)==0',None))
    constrainInvariantErrs()

    dvs = C.veccat([dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), dae.xDotVec()])
    ffcn = C.SXFunction([dvs],[sum([dae[n]**2 for n in ['aileron','elevator','y','z']])])
    gfcn = C.SXFunction([dvs],[g.getG()])
    ffcn.init()
    gfcn.init()

    guess = {'x':r0,'y':0,'z':0,
             'r':r0,'dr':0,
             'e11':0, 'e12':1, 'e13':0,
             'e21':0, 'e22':0, 'e23':1,
             'e31':1, 'e32':0, 'e33':0,
             'dx':0,'dy':0,'dz':0,
             'w1':0,'w2':omega0,'w3':0,
             'delta':0,'ddelta':omega0,
             'cos_delta':1,'sin_delta':0,
             'aileron':0,'elevator':0,
             'daileron':0,'delevator':0,
             'nu':300,'motor_torque':10,
             'ddr':0.0,'w0':0.0}
    dotGuess = {'x':0,'y':0,'z':0,'dx':0,'dy':0,'dz':0,
                'r':0,'dr':0,
                'e11':omega0,'e12':0,'e13':0,
                'e21':0,'e22':-omega0,'e23':0,
                'e31':0,'e32':0,'e33':0,
                'w1':0,'w2':0,'w3':0,
                'delta':omega0,'ddelta':0,
                'cos_delta':0,'sin_delta':omega0,
                'aileron':0,'elevator':0}

    guessVec = C.DMatrix([guess[n] for n in dae.xNames()+dae.zNames()+dae.uNames()+dae.pNames()]+
                         [dotGuess[n] for n in dae.xNames()])

    bounds = {'x':(0.01,r0+1),'y':(-1,0),'z':(-1,0),
             'r':(r0,r0),'dr':(0,0),
             'e11':(-2,2),'e12':(-2,2),'e13':(-2,2),
             'e21':(-2,2),'e22':(-2,2),'e23':(-2,2),
             'e31':(-2,2),'e32':(-2,2),'e33':(-2,2),
             'dx':(-50,50),'dy':(-50,50),'dz':(0,0),
             'w1':(-50,50),'w2':(-50,50),'w3':(-50,50),
             'delta':(0,0),'ddelta':(omega0,omega0),
             'cos_delta':(1,1),'sin_delta':(0,0),
             'aileron':(-0.1,0.1),'elevator':(-0.1,0.1),
             'daileron':(0,0),'delevator':(0,0),
             'nu':(0,3000),'motor_torque':(0,1000),
             'ddr':(0,0),'w0':(0,0)}
    dotBounds = {'x':(-50,50),'y':(-50,50),'z':(-50,50)
                 ,'dx':(0,0),'dy':(0,0),'dz':(0,0),
                 'r':(-1,1),'dr':(-1,1),
                 'e11':(-50,50),'e12':(-50,50),'e13':(-50,50),
                 'e21':(-50,50),'e22':(-50,50),'e23':(-50,50),
                 'e31':(-50,50),'e32':(-50,50),'e33':(-50,50),
                 'w1':(0,0),'w2':(0,0),'w3':(0,0),
                 'delta':(omega0-1,omega0+1),'ddelta':(0,0),
                 'cos_delta':(0,0),'sin_delta':(omega0,omega0),
                 'aileron':(-1,1),'elevator':(-1,1)}
    boundsVec = [bounds[n] for n in dae.xNames()+dae.zNames()+dae.uNames()+dae.pNames()]+ \
                [dotBounds[n] for n in dae.xNames()]
    

#    gfcn.setInput(guessVec)
#    gfcn.evaluate()
#    ret = gfcn.output()
#    for k,v in enumerate(ret):
#        if math.isnan(v):
#            print 'index ',k,' is nan: ',g._tags[k]
#    import sys; sys.exit()
    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")
    class MyCallback:
        def __init__(self):
            import rawekite.kiteproto as kiteproto
            import rawekite.kite_pb2 as kite_pb2
            self.kiteproto = kiteproto
            self.kite_pb2 = kite_pb2
            self.iter = 0
        def __call__(self,f,*args):
            x = C.DMatrix(f.input('x'))
            sol = {}
            for k,name in enumerate(dae.xNames()+dae.zNames()+dae.uNames()+dae.pNames()):
                sol[name] = x[k].at(0)
            lookup = lambda name: sol[name]
            kp = self.kiteproto.toKiteProto(lookup,
                                            lineAlpha=0.2)
            mc = self.kite_pb2.MultiCarousel()
            mc.horizon.extend([kp])
            mc.messages.append("iter: "+str(self.iter))
            self.iter += 1
            publisher.send_multipart(["multi-carousel", mc.SerializeToString()])

    
    solver = C.IpoptSolver(ffcn,gfcn)
    def addCallback():
        nd = len(boundsVec)
        nc = g.getLb().size()
        c = C.PyFunction( MyCallback(), C.nlpsolverOut(x=C.sp_dense(nd,1), f=C.sp_dense(1,1), lam_x=C.sp_dense(nd,1), lam_p = C.sp_dense(0,1), lam_g = C.sp_dense(nc,1), g = C.sp_dense(nc,1) ), [C.sp_dense(1,1)] )
        c.init()
        solver.setOption("iteration_callback", c)
#    addCallback()
    solver.setOption('max_iter',10000)
    solver.init()

    solver.setInput(g.getLb(),'lbg')
    solver.setInput(g.getUb(),'ubg')
    solver.setInput(guessVec,'x0')
    lb,ub = zip(*boundsVec)
    solver.setInput(C.DMatrix(lb), 'lbx')
    solver.setInput(C.DMatrix(ub), 'ubx')

    solver.solve()
    publisher.close()
    context.destroy()
    xOpt = solver.output('x')
    k = 0
    sol = {}
    for name in dae.xNames()+dae.zNames()+dae.uNames()+dae.pNames():
        sol[name] = xOpt[k].at(0)
        k += 1
#        print name+':\t',sol[name]
    dotSol = {}
    for name in dae.xNames():
        dotSol[name] = xOpt[k].at(0)
        k += 1
#        print 'DDT('+name+'):\t',dotSol[name]
    return sol

