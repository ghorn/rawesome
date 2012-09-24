import zmq
import time
import os

import casadi as C

import kite_pb2
import model
import joy

#tc0 = 2*389.970797939731

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

r = 1.2

ts = 0.02

def toProto(x,u):
    cs = kite_pb2.CarouselState()

    cs.kiteXyz.x = x.at(0)
    cs.kiteXyz.y = x.at(1)
    cs.kiteXyz.z = x.at(2)

    cs.kiteDcm.r11 = x.at(3)
    cs.kiteDcm.r12 = x.at(4)
    cs.kiteDcm.r13 = x.at(5)

    cs.kiteDcm.r21 = x.at(6)
    cs.kiteDcm.r22 = x.at(7)
    cs.kiteDcm.r23 = x.at(8)

    cs.kiteDcm.r31 = x.at(9)
    cs.kiteDcm.r32 = x.at(10)
    cs.kiteDcm.r33 = x.at(11)

    cs.delta = x.at(18)
    cs.ddelta = x.at(19)
    
    cs.tc = u.at(0)
    cs.u1 = u.at(1)
    cs.u2 = u.at(2)
    cs.wind_x = u.at(3)
    return cs

if __name__=='__main__':
    print "creating model"
    (dae, others) = model.model()
    dae.init()

    # compile model code
    print "generating model code"
    t0 = time.time()
    dae.generateCode("dae.c")
    print "took "+str(time.time()-t0)+" seconds to generate code"
    t0 = time.time()
    os.system("gcc -fPIC -O2 -shared dae.c -o dae.so")
    print "took "+str(time.time()-t0)+" seconds to compile code"
    dae_ext = C.ExternalFunction("./dae.so")
    dae_ext.init()
    
    print "creating integrator"
    f = C.IdasIntegrator(dae_ext)
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
        u1 = -axes[0]*0.02
        u2 =  axes[1]*0.05
        tc = 600*(1 - axes[6])
        wind_x = 5*(1-axes[7])

        u = C.DMatrix([tc,u1,u2,wind_x])
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
