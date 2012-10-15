import zmq
import time
import os
import math
import simutils

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

sim = simutils.Sim(ts=0.02, sloMoFactor=4, state0=simutils.SimState(pdOn = False, x = x0))

if __name__=='__main__':
    print "creating model"
    dae = model.model(zt)
    sxfun = dae.sxFun()
    sxfun.init()

    print "creating outputs function"
    fOutputs = dae.outputsFun()
    fOutputs.init()

    print "creating communicator"
    communicator = simutils.Communicator(fOutputs,dae._outputNames)

    print "creating integrator"
    f = C.IdasIntegrator(sxfun)
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
#        if len(js['buttonsDown']) > 0: print js['buttonsDown']

        # saving/loading states
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

        if 4 in js['buttonsDown']:
            print "saving file"
            sim.saveFile()
        if 2 in js['buttonsDown']:
            print "loading file"
            sim.loadFile()

        # play replay
        if 3 in js['buttonsDown']:
            sim.playReplay(communicator)
            sim.loadDefault()
        
        aileron = -js['axes'][0]*0.05
        elevator =  js['axes'][1]*0.2
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
        
        u = C.DMatrix([tc,aileron,elevator,ddr])
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
            except RuntimeError:
                sim.loadDefault()
                x,u,p = sim.getCurrentState()._log[-1]
                pass
            communicator.sendKite(sim,(x,u,p))
            
            deltaTime = (t0 + sim.tsSimStep*sim.sloMoFactor) - time.time()
            if deltaTime > 0:
                time.sleep(deltaTime)
    except KeyboardInterrupt:
        print "closing..."
        publisher.close()
        context.term()
        pass
