import zmq
import pickle

import casadi as C
from casadi import pi

from config import readConfig
import models
from newton.nmhe import Nmhe

if __name__=='__main__':
    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    print "reading config..."
#    conf = readConfig('stingray.ini','configspec.ini')
    conf = readConfig('config.ini','configspec.ini')
    conf['nk'] = 40
    
    print "creating model..."
    dae = models.crosswind(conf)#,extraParams=['endTime'])
    
    nk = 30

    # setup MHE
    mhe = Nmhe(dae,nk)

    # bounds
    mhe.bound('aileron',(-0.04,0.04))
    mhe.bound('elevator',(-0.1,0.3))

    mhe.bound('x',(-200,200))
    mhe.bound('y',(-200,200))
    if 'minAltitude' in conf:
        mhe.bound('z',(conf['minAltitude'],200))
    else:
        mhe.bound('z',(0.5,200))
    mhe.bound('r',(1,200))
    mhe.bound('dr',(-30,30))

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        mhe.bound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        mhe.bound(d,(-70,70))

    for w in ['w1','w2','w3']:
        mhe.bound(w,(-4*pi,4*pi))

#    mhe.bound('endTime',(0.5,10))
#    mhe.bound('endTime',(0.5,numLoops*7))
    mhe.bound('w0',(10,10))

    # boundary conditions
#    mhe.bound('y',(0,0),timestep=0)

    # guesses
#    mhe.guess('endTime',5.4)
    mhe.guess('w0',10)

    # constrain invariants
    def constrainInvariantErrs():
        dcm = mhe.lookup('dcm',timestep=0)
        err = C.mul(dcm.T,dcm)
        mhe.constrain( C.veccat([err[0,0] - 1, err[1,1]-1, err[2,2] - 1, err[0,1], err[0,2], err[1,2]]), '==', 0)
        mhe.constrain(mhe.lookup('c',timestep=0), '==', 0)
        mhe.constrain(mhe.lookup('cdot',timestep=0), '==', 0)
    constrainInvariantErrs()

    # initial guess
    print "loading trajectory..."
    f=open('data/crosswind_opt.dat','r')
    traj = pickle.load(f)
    f.close()
    
    for name in dae.xNames():
        for k in range(nk+1):
            mhe.guess(name,traj.lookup(name,timestep=k),timestep=k)
    for name in dae.pNames():
        mhe.guess(name,traj.lookup(name))

    # make objective
    obj = 0
    for name in ['x','y','z','dx','dy','dz','w1','w2','w3','aileron','elevator']:
        for k in range(nk+1):
            f = mhe.lookup(name,timestep=k) - traj.lookup(name,timestep=k)
            mhe.addGaussNewtonObjF(f)
#            obj += 0.5*f*f
#    for name in dae.pNames():
#        mhe.guess(name,traj.lookup(name))
#        nmhe.addGaussNewtonObjF(nmhe('x',timestep=k) - xTraj[k][0])
#        nmhe.addGaussNewtonObjF(nmhe('z',timestep=k) - xTraj[k][1])
#        nmhe.addGaussNewtonObjF(nmhe('dx',timestep=k) - xTraj[k][2])
#        nmhe.addGaussNewtonObjF(nmhe('dz',timestep=k) - xTraj[k][3])

#        obj += (nmhe('x',timestep=k) - xTraj[k][0])**2
#        obj += (nmhe('z',timestep=k) - xTraj[k][1])**2
#        obj += (nmhe('dx',timestep=k) - xTraj[k][2])**2
#        obj += (nmhe('dz',timestep=k) - xTraj[k][3])**2
        
#        obj += 1e-8*nmhe('x',timestep=k)**2
#        obj += 1e-8*nmhe('z',timestep=k)**2
#        obj += 1e-8*nmhe('dx',timestep=k)**2
#        obj += 1e-8*nmhe('dz',timestep=k)**2
#    obj += 1e-8*nmhe('m')**2
    mhe.setObj(obj)


#    nmhe.constrain(nmhe('dz',timestep=0)**2,'<=',1000)

    # get u
    uTraj = C.DMatrix([[traj.lookup(name,timestep=k) for k in range(nk)] for name in dae.uNames()]).T
    endTime = traj.tgrid[1,0,0] - traj.tgrid[0,0,0]
    mhe.makeSolver(endTime,traj=traj)
    mhe.runSolver(uTraj,traj)

#    # callback function
#    class MyCallback:
#        def __init__(self):
#            self.iter = 0 
#        def __call__(self,f,*args):
#            self.iter = self.iter + 1
#            xOpt = numpy.array(f.input(C.NLP_X_OPT))
#
#            traj = trajectory.Trajectory(ocp,xOpt)
#            
#            kiteProtos = []
#            for k in range(0,mhe.nk):
#                for nicpIdx in range(0,mhe.nicp):
#                    for j in [0]:
##                    for j in range(mhe.deg+1):
#                        kiteProtos.append( kiteproto.toKiteProto(C.DMatrix(traj.dvMap.xVec(k,nicpIdx=nicpIdx,degIdx=j)),
#                                                                 C.DMatrix(traj.dvMap.uVec(k)),
#                                                                 C.DMatrix(traj.dvMap.pVec()),
#                                                                 conf['kite']['zt'],
#                                                                 conf['carousel']['rArm'],
#                                                                 lineAlpha=0.2,
#                                                                 zeroDelta=True) )
#            mc = kite_pb2.MultiCarousel()
#            mc.css.extend(list(kiteProtos))
#
#            mc.messages.append("w0: "+str(traj.lookup('w0')))
#            mc.messages.append("iter: "+str(self.iter))
#            mc.messages.append("endTime: "+str(traj.lookup('endTime')))
#            mc.messages.append("average power: "+str(traj.lookup('quadrature energy',timestep=-1)/traj.lookup('endTime'))+" W")
#
#            # bounds feedback
##            lbx = mhe.solver.input(C.NLP_LBX)
##            ubx = mhe.solver.input(C.NLP_UBX)
##            violations = boundsFeedback(xOpt,lbx,ubx,mhe.bndtags,tolerance=1e-9)
##            for name in violations:
##                violmsg = "violation!: "+name+": "+str(violations[name])
##                mc.messages.append(violmsg)
#            
#            publisher.send_multipart(["multi-carousel", mc.SerializeToString()])


#    # Plot the results
#    def plotResults():
#        traj.subplot(['x','y','z'])
##        traj.subplot(['dx','dy','dz'])
#        traj.subplot([['aileron','elevator'],['daileron','delevator']],title='control surfaces')
#        traj.subplot(['r','dr','ddr'])
#        traj.subplot(['wind at altitude','dr'],title='')
#        traj.subplot(['c','cdot','cddot'],title="invariants")
#        traj.plot('airspeed',title='airspeed')
#        traj.subplot([['alpha(deg)','alphaTail(deg)'],['beta(deg)','betaTail(deg)']])
#        traj.subplot(['cL','cD','L/D'],title='')
#        traj.subplot([['winch power'], ['tether tension'],['accel (g)','accel without gravity (g)']])
#        traj.subplot([['ddx','ddy','ddz'],['accel','accel without gravity']])
#        traj.plot(["loyd's limit","loyd's limit (exact)","-(winch power)"])
#        traj.plot(["loyd's limit","-(winch power)"],title='')
#        traj.subplot([['daileronCost','delevatorCost','ddrCost'],['winch power']])
##        traj.subplot(['w1','w2','w3'])
#        traj.subplot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
#        traj.subplot(['energy','quadrature energy'])
#        traj.plot(['energy','quadrature energy'])
#        
#        plt.show()
#    plotResults()
