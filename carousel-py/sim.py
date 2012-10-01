import zmq
import time
import os
import math
import copy

import casadi as C

import kiteproto
import model
import joy

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
                ])

zt = -0.01

def setC0(x00):
    x = x00[0]
    y = x00[1]
    z = x00[2]
    e31 = x00[9]
    e32 = x00[10]
    e33 = x00[11]
    r = C.sqrt((x + zt*e31)**2 + (y + zt*e32)**2 + (z + zt*e33)**2)
    return C.veccat([x00,r,0])
    
x0 = setC0(x0)

class Communicator():
    def __init__(self, fOutputs, outputNames):
        self.context   = zmq.Context(1)
        self.publisher = self.context.socket(zmq.PUB)
        self.publisher.bind("tcp://*:5563")
        self.fOutputs = fOutputs
        self.outputNames = outputNames

    def sendKite(self,x,u,p,otherMessages=[]):
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
#        self.fOutputs.output()
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
            sim.sloMoFactor = 1 + 6*thumb
        else:
            sim.sloMoFactor = 1 + thumb*0.9
        return js

    def playReplay(self):
        print "replaying save #"+str(self.default)
        log = self._saves[self.default]._log
        loglen = float(len(log))
        for k,(x,u,p) in enumerate(log):
            t0 = time.time()
            self.handleInput() # for time dialation
            percent = "( %.1f %%)" % (100*k/loglen)
            message = "============ REPLAY #"+str(self.default)+" "+percent+" ==========="
            communicator.sendKite(x,u,p,[message])
            deltaTime = (t0 + sim.tsSimStep*sim.sloMoFactor) - time.time()
            if deltaTime > 0:
                time.sleep(deltaTime)



sim = Sim(ts=0.02, sloMoFactor=4, state0=SimState(pdOn = False, x = x0))

if __name__=='__main__':
    print "creating model"
    (dae, others, outputs) = model.model(zt)
    dae.init()

    print "creating outputs function"
    outputNames = outputs.keys()
    fOutputs = C.SXFunction([others['xVec'],C.veccat([others['uVec'],others['pVec']])],[outputs[n] for n in outputNames])
    fOutputs.setOption('name','fOutputs')
    fOutputs.init()
    communicator = Communicator(fOutputs,outputNames)

    print "creating integrator"
    f = C.IdasIntegrator(dae)
    f.setOption("reltol",1e-6)
    f.setOption("abstol",1e-8)
    f.setOption("t0",0)
    f.setOption("tf",sim.tsSimStep)
    f.setOption('name','integrator')
    f.setOption("linear_solver_creator",C.CSparse)
    f.setOption("linear_solver","user_defined")
    f.init()
    
    pd = {}
    pd['d'] = 0.1
    pd['p'] = 1e5
    
    def advanceState():
        js = sim.handleInput()
#        print [(k,v) for k,v in enumerate(js['axes'])]
#        print js['buttonsDown']

        # saving/loading
        fstButton=13
        for k in range(0,4):
            if k+fstButton+4 in js['buttonsDown'] and k in sim._saves:
                print "loading save #"+str(k+1)
                sim.load(k)
                sim.default = k
            if k+fstButton in js['buttonsDown']:
                if k==0:
                    print "can't override save0"
                else:
                    print "making save #"+str(k+1)
                    sim.save(k)
                    sim.default=k
        if 5 in js['buttonsDown']:
            k = sim.default
            print "loading save #"+str(k+1)
            sim.load(k)

        # play replay
        if 3 in js['buttonsDown']:
            sim.playReplay()
            sim.loadDefault()
        
        u1 = -js['axes'][0]*0.05
        u2 =  js['axes'][1]*0.2
        tc = 500*(1 - js['axes'][6])
        w0 = 10*(1-js['axes'][7])

        # torque PD controller
        if 9 in js['buttonsDown']:
            sim.currentState.pdOn = True
        if 10 in js['buttonsDown']:
            sim.currentState.pdOn = False

        x = sim.getCurrentState().x
        if sim.getCurrentState().pdOn:
            ddelta = x[19]
            delta = ((x[18].at(0)+math.pi) % (2*math.pi)) - math.pi # bound between (-pi,pi)
            tc = -pd['p']*(delta + pd['d']*ddelta)
#            print "PD: pgain: "+str(pd['p'])+" dgain: "+str(pd['d'])+" delta: "+str(delta)+" ddelta: "+str(ddelta)+" torque: "+str(tc)

        drRef = 5*js['axes'][10]

        r = x[20]
        dr = x[21]
        ddr = (drRef-dr)/sim.tsSimStep
#        print "r: %.2f\tdr: %.2f\tdrRef: %.2f\tdrErr: %.2f\tddr: %.2f" % (r,dr,drRef,drRef-dr,ddr)
        
        u = C.DMatrix([tc,u1,u2,ddr])
        p = C.DMatrix([w0])
        
        f.setInput(x,C.INTEGRATOR_X0)
        f.setInput(C.veccat([u,p]),C.INTEGRATOR_P)
        f.evaluate()

        xNext = C.DMatrix(f.output())
        return ((x, u, p), xNext)

    print "simulating..."
    try:
        while True:
            t0 = time.time()
            try:
                ((x,u,p), xNext) = advanceState()
                sim.currentState.log(x,u,p)
                sim.currentState.x = xNext
                if len(sim._saves[sim.default]._log)==0:
                    sim.save(sim.default)

                if len(sim.getCurrentState()._log)==0:
                    sim.log(x,u,p)
            except RuntimeError:
                sim.loadDefault()
                x,u,p = sim.getCurrentState()._log[-1]
                pass
            communicator.sendKite(x,u,p)
            
            deltaTime = (t0 + sim.tsSimStep*sim.sloMoFactor) - time.time()
            if deltaTime > 0:
                time.sleep(deltaTime)
    except KeyboardInterrupt:
        print "closing..."
        publisher.close()
        context.term()
        pass
