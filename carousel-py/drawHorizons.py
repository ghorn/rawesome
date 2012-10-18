import zmq
import kite_pb2
import time
import numpy

# zero mq setup
context   = zmq.Context(1)
publisher = context.socket(zmq.PUB)
publisher.bind("tcp://*:5563")

filename = "fullstate_mhe_mpc.dat"
f = open(filename,'r')
mhes = []
mpcs = []

nstates = 13
for line in f:
    x = map(float,line.split("\t"))
    x = [x[k:k+nstates] for k in range(0,len(x),nstates)]
    mhes.append( x[:11] )
    mpcs.append( x[11:] )
f.close()


def toKiteProto(x,alpha=1):
    cs = kite_pb2.CarouselState()

    cs.kiteXyz.x = x[0]
    cs.kiteXyz.y = x[1]
    cs.kiteXyz.z = x[2]

    cs.kiteDcm.r11 = x[3]
    cs.kiteDcm.r21 = x[4]
    cs.kiteDcm.r31 = x[5]

    cs.kiteDcm.r12 = x[6]
    cs.kiteDcm.r22 = x[7]
    cs.kiteDcm.r32 = x[8]

    cs.kiteDcm.r13 = x[9]
    cs.kiteDcm.r23 = x[10]
    cs.kiteDcm.r33 = x[11]

    cs.delta = x[12]

    cs.rArm = 1.08
    cs.zt = -0.01

    cs.transparency = alpha

    return cs

for k,(mhe,mpc) in enumerate(zip(mhes,mpcs)):

    n = len(mhe)
    
    mheProtos = []
    mpcProtos = []
    alphas = list(numpy.linspace(0.2,1,n))
    
    mheProtos = [toKiteProto(x,alpha=alpha) for alpha,x in zip(alphas,mhe[:n])]
    mpcProtos = [toKiteProto(x,alpha=alpha) for alpha,x in zip(reversed(alphas),mpc[:n])]
    
    mc = kite_pb2.MultiCarousel()
    mc.css.extend(list(mheProtos)+list(mpcProtos))

    mc.messages.append("number: "+str(k))
    publisher.send_multipart(["multi-carousel", mc.SerializeToString()])
    time.sleep(0.1)
