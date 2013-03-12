import zmq
from multiprocessing import Manager, Process
import time
import numpy
import casadi as C

from collocation import trajectory
import kiteproto
import kite_pb2

def normalCallback(traj,myiter,ocp,conf):
    kiteProtos = []
    for k in range(0,ocp.nk):
        for nicpIdx in range(0,ocp.nicp):
            for degIdx in [0]:
#            for degIdx in range(ocp.deg+1):
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
    

def startKiteTelemetry(ocp, conf,userCallback=normalCallback, printBoundViolation=False,printConstraintViolation=False):
    xOptQueue = Manager().Queue()
    KiteTelemetry(ocp, xOptQueue, conf, userCallback).start()

    class MyCallback:
        def __call__(self,f,*args):
            xOpt = numpy.array(f.input(C.NLP_X_OPT))

            if printBoundViolation:
                lbx = numpy.array(ocp.solver.input(C.NLP_LBX))
                ubx = numpy.array(ocp.solver.input(C.NLP_UBX))
                ocp._bounds.printBoundsFeedback(xOpt,lbx,ubx,reportThreshold=0)

            if printConstraintViolation:
                lbg = numpy.array(ocp.solver.input(C.NLP_LBG))
                ubg = numpy.array(ocp.solver.input(C.NLP_UBG))
                ocp._gfcn.setInput(xOpt,0)
                ocp._gfcn.evaluate()
                g = ocp._gfcn.output()
                s2 = ocp._constraints.printViolations(g,lbg,ubg,reportThreshold=0)

            xOptQueue.put(xOpt)

    return MyCallback()


class KiteTelemetry(Process):
    def __init__(self, ocp, xOptQueue, conf, callbackFun):
        Process.__init__(self)
        self.daemon = True

        self.ocp = ocp
        self.xOptQueue = xOptQueue
        self.conf = conf
        self.callbackFun = callbackFun

    def run(self):
        myiter = 0
        context   = zmq.Context(1)
        publisher = context.socket(zmq.PUB)
        publisher.bind("tcp://*:5563")

        while True:
            xOpt = self.xOptQueue.get()
            traj = trajectory.Trajectory(self.ocp,xOpt)
            mcStr = self.callbackFun(traj,myiter,self.ocp,self.conf)
            publisher.send_multipart(["multi-carousel", mcStr])
            myiter += 1
