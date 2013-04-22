import zmq
from multiprocessing import Manager, Process
import time
import numpy
import casadi as C

from rawe.collocation import trajectory
import kiteproto
import kite_pb2

def showAllPoints(traj,myiter,ocp,conf):
    return normalCallback(traj,myiter,ocp,conf,showAllPoints=True)

def normalCallback(traj,myiter,ocp,conf,showAllPoints=False):
    kiteProtos = []
    for k in range(0,ocp.nk):
        for nicpIdx in range(0,ocp.nicp):
            if showAllPoints:
                degIdxRange = range(ocp.deg+1)
            else:
                degIdxRange = [0]
            for degIdx in degIdxRange:
                lookup = lambda name: traj.lookup(name,timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)
                kiteProtos.append( kiteproto.toKiteProto(lookup,
                                                         lineAlpha=0.2) )
    t2 = time.time()
    mc = kite_pb2.MultiCarousel()
    mc.horizon.extend(list(kiteProtos))

    mc.messages.append("w0: "+str(traj.lookup('w0')))
    mc.messages.append("iter: "+str(myiter))
    mc.messages.append("endTime: "+str(traj.lookup('endTime')))
    mc.messages.append("average power: "+str(traj.lookup('quadrature energy',timestep=-1)/traj.lookup('endTime'))+" W")
    
    return mc.SerializeToString()
    

def startKiteTelemetry(ocp, conf, userCallback=normalCallback, printBoundViolation=False,printConstraintViolation=False, zeromqChannel='multi-carousel'):
    xOptQueue = Manager().Queue()
    KiteTelemetry(ocp, xOptQueue, conf, userCallback, zeromqChannel).start()

    class MyCallback:
        def __call__(self,f,*args):
            xOpt = numpy.array(f.input('x'))

            if printBoundViolation:
                lbx = numpy.array(ocp.solver.input('lbx'))
                ubx = numpy.array(ocp.solver.input('ubx'))
                ocp._bounds.printBoundsFeedback(xOpt,lbx,ubx,reportThreshold=0)

            if printConstraintViolation:
                lbg = numpy.array(ocp.solver.input('lbg'))
                ubg = numpy.array(ocp.solver.input('ubg'))
                ocp._gfcn.setInput(xOpt,0)
                ocp._gfcn.evaluate()
                g = ocp._gfcn.output()
                s2 = ocp._constraints.printViolations(g,lbg,ubg,reportThreshold=0)

            xOptQueue.put(xOpt)

    return MyCallback()


class KiteTelemetry(Process):
    def __init__(self, ocp, xOptQueue, conf, callbackFun, zeromqChannel):
        Process.__init__(self)
        self.daemon = True

        self.ocp = ocp
        self.xOptQueue = xOptQueue
        self.conf = conf
        self.callbackFun = callbackFun
        self.zeromqChannel = zeromqChannel

    def run(self):
        myiter = 0
        context   = zmq.Context(1)
        publisher = context.socket(zmq.PUB)
        publisher.bind("tcp://*:5563")

        while True:
            xOpt = self.xOptQueue.get()
            # empty the queue
            while not self.xOptQueue.empty():
                xOpt = self.xOptQueue.get()
            traj = trajectory.Trajectory(self.ocp,xOpt)
            mcStr = self.callbackFun(traj,myiter,self.ocp,self.conf)
            publisher.send_multipart([self.zeromqChannel, mcStr])
            myiter += 1
