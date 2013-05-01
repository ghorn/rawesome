import zmq
from multiprocessing import Manager, Process
import numpy

from rawe.collocation import trajectory

def startTelemetry(ocp, conf, callbacks=[],
                   printBoundViolation=False,printConstraintViolation=False,
                   url="tcp://*:5563"):
    xOptQueue = Manager().Queue()
    Sender(ocp, xOptQueue, conf, callbacks, url).start()

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


class Sender(Process):
    def __init__(self, ocp, xOptQueue, conf, callbackFunsWithChannels, url):
        Process.__init__(self)
        self.daemon = True

        self.ocp = ocp
        self.xOptQueue = xOptQueue
        self.conf = conf
        self.url = url
        msg = "callbacks must be a list of (callback (function), channel (string)) tuples"
        if not isinstance(callbackFunsWithChannels, list):
            callbackFunsWithChannels = [callbackFunsWithChannels]
        for cfc in callbackFunsWithChannels:
            assert isinstance(cfc, tuple), msg
            assert isinstance(cfc[1], basestring), msg
        self.callbackFunsWithChannels = callbackFunsWithChannels

    def run(self):
        myiter = 0
        context   = zmq.Context(1)
        publisher = context.socket(zmq.PUB)
        publisher.bind(self.url)

        while True:
            xOpt = self.xOptQueue.get()
            # empty the queue
            while not self.xOptQueue.empty():
                xOpt = self.xOptQueue.get()
                myiter += 1
            traj = trajectory.Trajectory(self.ocp,xOpt)
            for callbackFun, zeromqChannel in self.callbackFunsWithChannels:
                mcStr = callbackFun(traj,myiter,self.ocp,self.conf)
                publisher.send_multipart([zeromqChannel, mcStr])
