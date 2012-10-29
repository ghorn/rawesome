import matplotlib.pyplot as plt
import pickle
import numpy as np
import sys

import zmq
import kite_pb2
import kiteproto

def toKiteProto(x,k,zt=0,transparency=1,dz=0):
    cs = kite_pb2.CarouselState()

    cs.kiteXyz.x = x['x'][k]
    cs.kiteXyz.y = x['y'][k]
    cs.kiteXyz.z = x['z'][k]+dz

    cs.kiteDcm.r11 = x['e11'][k]
    cs.kiteDcm.r12 = x['e12'][k]
    cs.kiteDcm.r13 = x['e13'][k]

    cs.kiteDcm.r21 = x['e21'][k]
    cs.kiteDcm.r22 = x['e22'][k]
    cs.kiteDcm.r23 = x['e23'][k]

    cs.kiteDcm.r31 = x['e31'][k]
    cs.kiteDcm.r32 = x['e32'][k]
    cs.kiteDcm.r33 = x['e33'][k]

    if 'delta' in x:
        cs.delta = x['delta'][k]
    else:
        cs.delta = 0

    cs.rArm = 0
    cs.zt = zt

    cs.transparency = transparency

    return cs

context   = zmq.Context(1)
publisher = context.socket(zmq.PUB)
publisher.bind("tcp://*:5563")

# load saved trajectory
f=open("crosswind_opt.dat",'r')
opt = pickle.load(f)
f.close()

names     = opt['vardict'].keys()+['yaw','pitch','roll']

cosOrder = range(0,14)
sinOrder = range(1,14)

def fitFourier(ts_,xs):
    ts = ts_/ts_[-1]*2*np.pi
    M = []
    for k in range(ts.size):
        row = []
        t = ts[k]
        x = xs[k]
        for so in cosOrder:
            row.append(np.cos(so*t))
        for so in sinOrder:
            row.append(np.sin(so*t))
        M.append( np.array(row) )
    M = np.array( M )

    fitcoeffs = np.linalg.lstsq(M,xs)[0]
    return (fitcoeffs, np.dot(M,fitcoeffs))

def dcm2Euler(k):
    r11 = opt['vardict']['e11'][k]
    r12 = opt['vardict']['e12'][k]
    mr13 = -opt['vardict']['e13'][k]
    #mr13 -- nan protect
    #  | mr13' >  1 =  1
    #  | mr13' < -1 = -1
    #  | otherwise = mr13'
    r23 = opt['vardict']['e23'][k]
    r33 = opt['vardict']['e33'][k]
          
    yaw   = np.arctan2(r12,r11)
    pitch = np.arcsin(mr13)
    roll  = np.arctan2(r23,r33)
    return (yaw,pitch,roll)

# convert dcms to eulers
yaw = []
pitch = []
roll = []
for k in range(0,opt['tgrid'].size):
    (y,p,r) = dcm2Euler(k)
    yaw.append(float(y))
    pitch.append(float(p))
    roll.append(float(r))
opt['vardict']['yaw']   = np.array(yaw)
opt['vardict']['pitch'] = np.array(pitch)
opt['vardict']['roll']  = np.array(roll)

# fit everything
fitcoeffs = {}
fit = {}

sys.stdout.write("fitting: ")
sys.stdout.flush()
for name in names:
    # don't fit parameters
    if not isinstance(opt['vardict'][name], float):
        sys.stdout.write(name+" ")
        sys.stdout.flush()
        (fitcoeffs_,fit_) = fitFourier(opt['tgrid'], opt['vardict'][name])
        fitcoeffs[name] = fitcoeffs_
        fit[name] = fit_
sys.stdout.write('\n')
sys.stdout.flush()

# send kite protos
kiteProtos = []
for k in range(0,opt['tgrid'].size):
    kiteProtos.append( toKiteProto(opt['vardict'],k,zt=0,transparency=0.2) )
    kiteProtos.append( toKiteProto(fit,k,zt=0, transparency=1,dz=1) )
mc = kite_pb2.MultiCarousel()
mc.css.extend(list(kiteProtos))
publisher.send_multipart(["multi-carousel", mc.SerializeToString()])

def plot(plotnames):
    plt.figure()
    plt.clf()
    legend = []
    for name in plotnames:
        legend.append(name)
        plt.plot(opt['tgrid'],opt['vardict'][name])
        legend.append(name+" fit")
        plt.plot(opt['tgrid'],fit[name])
    plt.title("fourier fits")
    plt.xlabel('time')
    plt.legend(legend)
    plt.grid()
#plot(['x','y','z'])
#plot(['dx','dy','dz'])
plot(['w1','w2','w3'])
plot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
#plot(['yaw','pitch','roll'])
#plot(['r'])
#plot(['dr'])
plt.show()

fitcoeffs['sinOrder'] = sinOrder
fitcoeffs['cosOrder'] = cosOrder

print "saving fourier coeffs"
f=open("crosswind_opt_fourier.dat",'w')
pickle.dump(opt,f)
f.close()
