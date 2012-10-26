import zmq
import time
import os
import math
import copy
import pickle

import casadi as C

import kiteproto
import joy

class Communicator():
    def __init__(self, fOutputs, outputNames):
        self.context   = zmq.Context(1)
        self.publisher = self.context.socket(zmq.PUB)
        self.publisher.bind("tcp://*:5563")
        self.fOutputs = fOutputs
        self.outputNames = outputNames

    def sendKite(self,sim,(x,z,u,p),zt,rArm,w0,otherMessages=[]):
        assert(isinstance(otherMessages,list))
        pb = kiteproto.toKiteProto(x,u,p,zt,rArm,w0=w0)
        pb.messages.append("sloMoFactor: "+str(sim.sloMoFactor))
        pb.messages.append("-------------------------")

        self.fOutputs.setInput(x,0)
        self.fOutputs.setInput(z,1)
        self.fOutputs.setInput(u,2)
        self.fOutputs.setInput(p,3)
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
        assert( x is not None )
        
        self.pdOn=pdOn
        self.x = x
        self._log = []

    def log(self,x,u,p):
        self._log.append((x,u,p))


class Sim():
    def __init__(self,ts=None, sloMoFactor=None, state0=None):
        assert(ts is not None)
        assert(sloMoFactor is not None)
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

    def saveFile(self,filename='defaultSave.dat'):
        saves = copy.deepcopy(self._saves)
        for k in saves.keys():
            if not k == 'default':
                saves[k].x = list(saves[k].x)
                saves[k]._log = [(list(x),list(u),list(p)) for (x,u,p) in saves[k]._log]
        saves['default'] = self.default

        f=open(filename,'w')
        pickle.dump(saves,f)
        f.close()

    def loadFile(self,filename='defaultSave.dat'):
        f=open(filename,'r')
        saves = pickle.load(f)
        f.close()
        
        for k in saves.keys():
            if isinstance(k,int):
                saves[k].x = C.DMatrix(saves[k].x)
                saves[k]._log = [(C.DMatrix(x),C.DMatrix(u),C.DMatrix(p)) for (x,u,p) in saves[k]._log]
        self.default = saves['default']
        self._saves = saves
        self.load(self.default)

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
            communicator.sendKite(self,(x,C.DMatrix([0,0,0,0,0,0,0,0]),u,p),0.0,1,0,[message])
            deltaTime = (t0 + self.tsSimStep*self.sloMoFactor) - time.time()
            if deltaTime > 0:
                time.sleep(deltaTime)
