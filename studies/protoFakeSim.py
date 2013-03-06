import pickle
from scipy.interpolate import PiecewisePolynomial
import zmq

import rawe

def loadPps(filename):
    print "loading ",filename
    f=open(filename,'r')
    traj = pickle.load(f)
    f.close()

    assert isinstance(traj,rawe.collocation.trajectory.Trajectory), "the file \""+filename+"\" doean't have a pickled Trajectory"

#    h = (traj.tgrid[-1,0,0] - traj.tgrid[0,0,0])/float(traj.dvMap._nk*traj.dvMap._nicp)
#    h *= traj.dvMap._nk*traj.dvMap._nicp/float(self.nk*self.nicp)
#    h *= numLoops
        
    pps = {}
#    missing = []
    ############# make piecewise polynomials ###########
    # differential states
    for name in traj.dvMap._xNames:
        # make piecewise poly
        pps[name] = None
        for timestepIdx in range(traj.dvMap._nk):
            for nicpIdx in range(traj.dvMap._nicp):
                ts = []
                ys = []
                for degIdx in range(traj.dvMap._deg+1):
                    ts.append(traj.tgrid[timestepIdx,nicpIdx,degIdx])
                    ys.append([traj.dvMap.lookup(name,timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx)])
                if pps[name] is None:
                    pps[name] = PiecewisePolynomial(ts,ys)
                else:
                    pps[name].extend(ts,ys)
        pps[name].extend([traj.tgrid[-1,0,0]],[[traj.dvMap.lookup(name,timestep=-1,nicpIdx=0,degIdx=0)]])

    # algebraic variables
    for name in traj.dvMap._zNames:
        # make piecewise poly
        pps[name] = None
        for timestepIdx in range(traj.dvMap._nk):
            for nicpIdx in range(traj.dvMap._nicp):
                ts = []
                ys = []
                for degIdx in range(1,traj.dvMap._deg+1):
                    ts.append(traj.tgrid[timestepIdx,nicpIdx,degIdx])
                    ys.append([traj.dvMap.lookup(name,timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx)])
                if pps[name] is None:
                    pps[name] = PiecewisePolynomial(ts,ys)
                else:
                    pps[name].extend(ts,ys)

    # controls
    for name in traj.dvMap._uNames:
        # make piecewise poly
        ts = []
        ys = []
        for timestepIdx in range(traj.dvMap._nk):
            ts.append(traj.tgrid[timestepIdx,0,0])
            ys.append([traj.dvMap.lookup(name,timestep=timestepIdx)])
        pps[name] = PiecewisePolynomial(ts,ys)


    return (traj,pps)


#    ############# interpolate ###########
#    # interpolate differential states
#    for name in self.dae.xNames():
#        if name not in pps:
#            missing.append(name)
#            continue
#        # evaluate piecewise poly to set initial guess
#        t0 = 0.0
#        for timestepIdx in range(self.nk):
#            for nicpIdx in range(self.nicp):
#                for degIdx in range(self.deg+1):
#                    time = t0 + h*self.lagrangePoly.tau_root[degIdx]
#                    if time > traj.tgrid[-1,0,0]:
#                        time -= traj.tgrid[-1,0,0]
#                    self.guess(name,pps[name](time),timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx,force=force,quiet=quiet)
#                t0 += h
#                if t0 > traj.tgrid[-1,0,0]:
#                    t0 -= traj.tgrid[-1,0,0]
#        self.guess(name,pps[name](t0),timestep=-1,nicpIdx=0,degIdx=0,force=force,quiet=quiet)
#
#
#    # interpolate algebraic variables
#    for name in self.dae.zNames():
#        if name not in pps:
#            missing.append(name)
#            continue
#        # evaluate piecewise poly to set initial guess
#        t0 = 0.0
#        for timestepIdx in range(self.nk):
#            for nicpIdx in range(self.nicp):
#                for degIdx in range(1,self.deg+1):
#                    time = t0 + h*self.lagrangePoly.tau_root[degIdx]
#                    if time > traj.tgrid[-1,0,0]:
#                        time -= traj.tgrid[-1,0,0]
#                    self.guess(name,pps[name](time),timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx,force=force,quiet=quiet)
#                t0 += h
#
#    # interpolate controls
#    for name in self.dae.uNames():
#        if name not in pps:
#            missing.append(name)
#            continue
#        # evaluate piecewise poly to set initial guess
#        t0 = 0.0
#        for timestepIdx in range(self.nk):
#            self.guess(name,pps[name](t0),timestep=timestepIdx,force=force,quiet=quiet)
#            t0 += h
#
#    # set parameters
#    for name in self.dae.pNames():
#        if name not in traj.dvMap._pNames:
#            missing.append(name)
#            continue
#        if name=='endTime':
#            self.guess(name,traj.dvMap.lookup(name)*numLoops,force=force,quiet=quiet)
#        else:
#            self.guess(name,traj.dvMap.lookup(name),force=force,quiet=quiet)
#
#    msg = "finished interpolating initial guess"
#    if len(missing) > 0:
#        msg += ", couldn't find fields: "+str(missing)
#    else:
#        msg += ", all fields found"
#    print msg

#filename = 'data/crosswind_opt_4_loops.dat'
filename = 'data/crosswind_opt.dat'
(traj,pps) = loadPps(filename)
dt = 1.0/30.0

# make a list of protos
print 'building list of protos'
protos = []
time = 0.0
while time < 100:
    t = time
    while t>traj.tgrid[-1,0,0]:
        t -= traj.tgrid[-1,0,0]
    lookup = lambda name: pps[name](t)
    kp = rawe.kiteproto.toKiteProto(lookup,lineAlpha=0.2)

    mc = rawe.kite_pb2.MultiCarousel()
    mc.css.extend([kp])
    mc.messages.append('time: %.3f' % time)

    protos.append( mc.SerializeToString() )
    time += dt

# send all the messages
context   = zmq.Context(1)
publisher = context.socket(zmq.PUB)
publisher.bind("tcp://*:5563")

import time
print 'sending protos'
k = 0
tNext = time.time() + dt
for mcStr in protos:
    publisher.send_multipart(["multi-carousel", mcStr])
    print "wooo ",k
    k += 1
    toSleep = tNext - time.time()
    if toSleep > 0:
        time.sleep(toSleep)
    tNext += dt
