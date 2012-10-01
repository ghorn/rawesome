import zmq
import time
import os
import math
import copy

import casadi as C

import kiteproto
import model
import joy

class Communicator():
    def __init__(self, fOutputs, outputNames):
        self.context   = zmq.Context(1)
        self.publisher = self.context.socket(zmq.PUB)
        self.publisher.bind("tcp://*:5563")
        self.fOutputs = fOutputs
        self.outputNames = outputNames

    def sendKite(self,sim,(x,u,p),otherMessages=[]):
        assert(isinstance(otherMessages,list))
        pb = kiteproto.toKiteProto(x,u,p)
        pb.messages.append("sloMoFactor: "+str(sim.sloMoFactor))
        pb.messages.append("torque: "+str(u.at(0)))
        pb.messages.append("u1: %.2f"%(u.at(1)*180/math.pi)+" deg")
        pb.messages.append("u2: %.2f"%(u.at(2)*180/math.pi)+" deg")
        pb.messages.append("r:  %.3f"%x.at(20))
        pb.messages.append("dr: %.3f"%x.at(21))
        pb.messages.append("RPM: "+str(x.at(19)*60/(2*math.pi)))
        pb.messages.append("-------------------------")

        self.fOutputs.setInput(x,0)
        self.fOutputs.setInput(C.veccat([u,p]),1)
        self.fOutputs.evaluate()
        for k,n in enumerate(self.outputNames):
            pb.messages.append(n+": "+str(self.fOutputs.output(k)))
        if len(otherMessages)>0:
            pb.messages.append("-------------------------")
            for om in otherMessages:
                assert(isinstance(om,str))
                pb.messages.append(om)

        self.publisher.send_multipart(["carousel", pb.SerializeToString()])
        

class SimState():
    def __init__(self,pdOn=None, x=None):
        assert( isinstance(pdOn, bool) )
        assert( x != None )
        
        self.pdOn=pdOn
        self.x = x
        self._log = []

    def log(self,x,u,p):
        self._log.append((x,u,p))


class Sim():
    def __init__(self,ts=None, sloMoFactor=None, state0=None):
        assert(ts!=None)
        assert(sloMoFactor!=None)
        assert(isinstance(state0, SimState))

        self.ts = ts
        self.sloMoFactor=sloMoFactor
        self.tsSimStep=self.ts/self.sloMoFactor
        
        self.currentState = state0
        self._saves = {}
        self.default = 0
        self.save(self.default)
        self.js = joy.Joy()

    def getCurrentState(self):
        return self.currentState

    def save(self,k):
        assert(isinstance(k,int))
        self._saves[k] = copy.deepcopy(self.getCurrentState())
#        print "log: "+str(self._saves[k]._log)
#        print "pdOn: "+str(self._saves[k].pdOn)
#        print "x: "+str(self._saves[k].x)

    def load(self,k):
        assert(isinstance(k,int))
        if k in self._saves:
            self.currentState = copy.deepcopy(self._saves[k])
#            print "log: "+str(self.currentState._log)
#            print "pdOn: "+str(self.currentState.pdOn)
#            print "x: "+str(self.currentState.x)

        else:
            print "can't load #"+str(k)+" because that save doesn't exist"

    def loadDefault(self):
        self.load(self.default)

    def handleInput(self):
        js = self.js.getAll()

        # time dialation
        thumb = js['axes'][8]
        if thumb>=0:
            self.sloMoFactor = 1 + 6*thumb
        else:
            self.sloMoFactor = 1 + thumb*0.9
        return js

    def playReplay(self, communicator):
        print "replaying save #"+str(self.default)
        log = self._saves[self.default]._log
        loglen = float(len(log))
        for k,(x,u,p) in enumerate(log):
            t0 = time.time()
            self.handleInput() # for time dialation
            percent = "( %.1f %%)" % (100*k/loglen)
            message = "============ REPLAY #"+str(self.default)+" "+percent+" ==========="
            communicator.sendKite(self,(x,u,p),[message])
            deltaTime = (t0 + self.tsSimStep*self.sloMoFactor) - time.time()
            if deltaTime > 0:
                time.sleep(deltaTime)
