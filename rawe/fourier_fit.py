import matplotlib.pyplot as plt
import pickle
import numpy as np
import sys
import copy

import zmq
import kite_pb2
import kiteproto

class FourierFit():
    def __init__(self,name,polyOrder_,cosOrder_,sinOrder_,ts_, xs_):
        self.name = name
        self.polyOrder = polyOrder_
        self.cosOrder  = cosOrder_
        self.sinOrder  = sinOrder_
        
        self.ts = copy.deepcopy(ts_)
        self.timeScaling = float(2*np.pi/self.ts[-1])
        self.scaledTs = self.ts*self.timeScaling
        self.xs = copy.deepcopy(xs_)
        self.M = []

        for k in range(self.ts.size):
            row = []
            t = self.scaledTs[k]
            x = self.xs[k]
            for po in self.polyOrder:
                row.append(t**po)
            for co in self.cosOrder:
                row.append(np.cos(co*t))
            for so in self.sinOrder:
                row.append(np.sin(so*t))
            self.M.append( np.array(row) )

        self.M = np.array( self.M )
        self.fitcoeffs = np.linalg.lstsq(self.M,self.xs)[0]
        self.Mx = np.dot(self.M,self.fitcoeffs)

        self.polyCoeffs = self.fitcoeffs[0:len(self.polyOrder)]
        self.cosCoeffs  = self.fitcoeffs[len(self.polyOrder):len(self.polyOrder+self.cosOrder)]
        self.sinCoeffs  = self.fitcoeffs[len(self.polyOrder+self.cosOrder):len(self.polyOrder+self.cosOrder+self.sinOrder)]
            
    def evaluate(self,t_):
        t = t_*self.timeScaling
        val = 0
        for k,po in enumerate(self.polyOrder):
            val += self.polyCoeffs[k]*t**po
        for k,so in enumerate(self.cosOrder):
            val += self.cosCoeffs[k]*np.cos(so*t)
        for k,so in enumerate(self.sinOrder):
            val += self.sinCoeffs[k]*np.sin(so*t)
        return val
        
if __name__=='__main__':
    filename = "data/crosswind_opt"
#    filename = "data/carousel_opt"
    
    def toKiteProto(fits,k,zt=0,kiteAlpha=1,lineAlpha=1, dz=0):
        cs = kite_pb2.CarouselState()
    
        cs.kiteXyz.x = fits['x'].Mx[k]
        cs.kiteXyz.y = fits['y'].Mx[k]
        cs.kiteXyz.z = fits['z'].Mx[k]+dz
    
        cs.kiteDcm.r11 = fits['e11'].Mx[k]
        cs.kiteDcm.r12 = fits['e12'].Mx[k]
        cs.kiteDcm.r13 = fits['e13'].Mx[k]
    
        cs.kiteDcm.r21 = fits['e21'].Mx[k]
        cs.kiteDcm.r22 = fits['e22'].Mx[k]
        cs.kiteDcm.r23 = fits['e23'].Mx[k]
    
        cs.kiteDcm.r31 = fits['e31'].Mx[k]
        cs.kiteDcm.r32 = fits['e32'].Mx[k]
        cs.kiteDcm.r33 = fits['e33'].Mx[k]
    
        if 'delta' in fits:
            cs.delta = fits['delta'].Mx[k]
        else:
            cs.delta = 0
    
        cs.rArm = 0
        cs.zt = zt
    
        cs.kiteTransparency = kiteAlpha
        cs.lineTransparency = lineAlpha
    
        return cs
    
    def npToKiteProto(x,k,zt=0,kiteAlpha=1,lineAlpha=1, dz=0):
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
    
        cs.kiteTransparency = kiteAlpha
        cs.lineTransparency = lineAlpha
    
        return cs
    
    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")
    
    # load saved trajectory
    loadfile = filename+".dat"
    print "loading fourier coeffs: "+loadfile
    f=open(loadfile,'r')
    opt = pickle.load(f)
    f.close()
    
    names = opt['vardict'].keys()+['yaw','pitch','roll']
    
    polyOrder = [0,1]
    cosOrder = range(1,14)
    sinOrder = range(1,14)
    
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
    fits = {}
    
    sys.stdout.write("fitting: ")
    sys.stdout.flush()
    for name in names:
        # don't fit parameters
        if not isinstance(opt['vardict'][name], float):
            sys.stdout.write(name+" ")
            sys.stdout.flush()
            ts = opt['tgrid']
            xs = opt['vardict'][name]
            fits[name] = FourierFit(name, polyOrder, cosOrder, sinOrder, ts, xs)
    sys.stdout.write('\n')
    sys.stdout.flush()
    
    # send kite protos
    kiteProtos = []
    for k in range(0,opt['tgrid'].size):
        kiteProtos.append( npToKiteProto(opt['vardict'],k,zt=0,kiteAlpha=0.2,lineAlpha=0.2) )
        kiteProtos.append( toKiteProto(fits,k,zt=0, dz=1) )
    mc = kite_pb2.MultiCarousel()
    mc.css.extend(list(kiteProtos))
    publisher.send_multipart(["multi-carousel", mc.SerializeToString()])
    
    def plot(plotnames):
        plt.figure()
        plt.clf()
        legend = []
        for name in plotnames:
            legend.append(name+" fit")
    
            t0 = opt['tgrid'][0]
            tf = opt['tgrid'][-1]
            dt = tf-t0
    
            t0 -= 0.2*dt
            tf += 0.2*dt
            ts = np.linspace(t0,tf,500)
            plt.plot(ts,[fits[name].evaluate(t) for t in ts])
    
            legend.append(name)
            plt.plot(opt['tgrid'],opt['vardict'][name],'--')
            
        plt.title("fourier fits")
        plt.xlabel('time')
        plt.legend(legend)
        plt.grid()
    
    savefile = filename+"_fourier.dat"
    print "saving fourier coeffs: "+savefile
    f=open(savefile,'w')
    pickle.dump(fits,f)
    f.close()
    
    def mkPlots():
        plot(['x','y','z'])
        plot(['dx','dy','dz'])
        plot(['w1','w2','w3'])
        plot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        plot(['yaw','pitch','roll'])
        plot(['r'])
        plot(['dr'])
        if 'delta' in opt['vardict']:
            plot(['delta'])
        if 'ddelta' in opt['vardict']:
            plot(['ddelta'])
        plt.show()
    mkPlots()
