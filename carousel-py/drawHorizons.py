import zmq
import kite_pb2
import time
import numpy

# zero mq setup
context   = zmq.Context(1)
publisher = context.socket(zmq.PUB)
publisher.bind("tcp://*:5563")

#filename = "fullstate_mhe_mpc.dat"
#filename = "fullstate_mhe_mpc_40rpm.dat"
#filename = "fullstate_mhe_mpc_45rpm_001.dat"
filename = "fullstate_mhe_mpc_45rpm_002.dat"
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

N0 = 85
NF = 130
k = N0
for (mhe,mpc) in zip(mhes[N0:NF],mpcs[N0:NF]):

    n = len(mhe)
#    n -= 5
    
    mheProtos = []
    mpcProtos = []
    alphas = list(numpy.linspace(0.1,0.4,n))
    alphas[-1] = 1

#    for alpha,x in zip(alphas,mhe[-n:]):
#        mheProtos.append( toKiteProto(x,alpha=alpha) )
#        print "mhe delta: "+str(x[12])
#        print alpha

    for alpha,x in zip(reversed(alphas),mpc[:n]):
        mpcProtos.append( toKiteProto(x,alpha=alpha) )
#        print "mpc delta: "+str(x[12])
#        print alpha
    
    mc = kite_pb2.MultiCarousel()
    mc.css.extend(list(mheProtos)+list(mpcProtos))

    mc.messages.append("number: "+str(k+1)+"/"+str(len(mhes)))
    publisher.send_multipart(["multi-carousel", mc.SerializeToString()])
    time.sleep(1)

    k += 1
