import zmq
import time
import os

import casadi as C

import kite_pb2
import pendulum_model
import joy

ts = 0.02
r = 0.3
x0 = [r,0,0,0]

def toProto(x,u):
    cs = kite_pb2.CarouselState()

    cs.kiteXyz.x = x.at(0)
    cs.kiteXyz.y = x.at(1)
    cs.kiteXyz.z = 0

    cs.kiteDcm.r11 = 1
    cs.kiteDcm.r12 = 0
    cs.kiteDcm.r13 = 0

    cs.kiteDcm.r21 = 0
    cs.kiteDcm.r22 = 1
    cs.kiteDcm.r23 = 0

    cs.kiteDcm.r31 = 0
    cs.kiteDcm.r32 = 0
    cs.kiteDcm.r33 = 1

    cs.delta = 0
    cs.ddelta = 0

    cs.tc = u.at(0)
    cs.u1 = 0
    cs.u2 = u.at(1)
    cs.wind_x = 0
    return cs

if __name__=='__main__':
    print "creating model"
    (dae, others) = pendulum_model.pendulum_model()
    dae.init()

#    # compile model code
#    print "generating model code"
#    t0 = time.time()
#    dae.generateCode("dae.c")
#    print "took "+str(time.time()-t0)+" seconds to generate code"
#    t0 = time.time()
#    os.system("gcc -fPIC -O2 -shared dae.c -o dae.so")
#    print "took "+str(time.time()-t0)+" seconds to compile code"
#    dae_ext = C.ExternalFunction("./dae.so")
#    dae_ext.init()
#    dae = dae_ext
    
    print "creating integrator"
    f = C.IdasIntegrator(dae)
    f.setOption("reltol",1e-6)
    f.setOption("abstol",1e-8)
    f.setOption("t0",0)
    f.setOption("tf",ts)
    f.init()
    
    js = joy.Joy()

    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    def advanceState(x):
        axes = js.getAxes()
        torque = axes[0]
        u = C.DMatrix([torque,0.4])
        
        f.setInput(x,C.INTEGRATOR_X0)
        f.setInput(u,C.INTEGRATOR_P)
        f.evaluate()
        return (C.DMatrix(f.output()), u)

    x = x0
    print "simulating..."
    try:
        while True:
            t0 = time.time()
            (x,u) = advanceState(x)
            p = toProto(x,u)
            publisher.send_multipart(["carousel", p.SerializeToString()])
            
            deltaTime = (t0 + ts) - time.time()
            if deltaTime > 0:
                time.sleep(deltaTime)

    except KeyboardInterrupt:
        print "closing..."
        publisher.close()
        context.term()
        pass
