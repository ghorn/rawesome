# Copyright 2012-2013 Greg Horn
#
# This file is part of rawesome.
#
# rawesome is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rawesome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib.pyplot as plt
import pickle
import numpy as np
import sys
import copy
import zmq

import casadi as C
import rawe

class TrajFit():
    def __init__(self,orderMap,traj):
        ts = [traj.tgrid[k,j,i] for k in range(traj.nk) for j in range(traj.nicp) for i in range(traj.deg+1)]
        ts.append(traj.tgrid[traj.nk,0,0])
        self.ts = np.array(ts)

        self.fits = {}

        def getAllX(name):
            xs = [traj.lookup(name,timestep=k,nicpIdx=j,degIdx=i) \
                  for k in range(traj.nk) for j in range(traj.nicp) for i in range(traj.deg+1)]
            xs.append(traj.lookup(name,timestep=traj.nk,nicpIdx=0,degIdx=0))
            xs = np.array(xs)
            return xs

        sys.stdout.write("fitting: ")
        sys.stdout.flush()
        for name in traj.dvMap._xNames:
            sys.stdout.write(name+" ")
            sys.stdout.flush()

            xs = getAllX(name)

            polyOrder = orderMap[name]['poly']
            sinOrder = orderMap[name]['sin']
            cosOrder = orderMap[name]['cos']
            self.fits[name] = FourierFit(name, polyOrder, cosOrder, sinOrder, self.ts, xs)

        sys.stdout.write('\n')
        sys.stdout.flush()

        # make kite protobufs
        t0 = self.ts[0]
        tf = self.ts[-1]
        dt = tf-t0
        ts_ = np.linspace(t0,tf,100)
        kiteProtos = []
        for t in ts_:
            cs = rawe.kite_pb2.CarouselState()

            cs.kiteXyz.x = self.fits['x'].evaluate(t)
            cs.kiteXyz.y = self.fits['y'].evaluate(t)
            cs.kiteXyz.z = self.fits['z'].evaluate(t)

            cs.kiteDcm.r11 = self.fits['e11'].evaluate(t)
            cs.kiteDcm.r12 = self.fits['e12'].evaluate(t)
            cs.kiteDcm.r13 = self.fits['e13'].evaluate(t)

            cs.kiteDcm.r21 = self.fits['e21'].evaluate(t)
            cs.kiteDcm.r22 = self.fits['e22'].evaluate(t)
            cs.kiteDcm.r23 = self.fits['e23'].evaluate(t)

            cs.kiteDcm.r31 = self.fits['e31'].evaluate(t)
            cs.kiteDcm.r32 = self.fits['e32'].evaluate(t)
            cs.kiteDcm.r33 = self.fits['e33'].evaluate(t)

            if 'delta' in self.fits:
                cs.delta = self.fits['delta'].evaluate(t)
            else:
                cs.delta = 0

            cs.rArm = 0
            cs.zt = 0

            cs.kiteTransparency = 0.2
            cs.lineTransparency = 0.2

            kiteProtos.append(cs)
        mc = rawe.kite_pb2.MultiCarousel()
        mc.horizon.extend(list(kiteProtos))
        self.multiCarousel = mc.SerializeToString()
        #publisher.send_multipart(["multi-carousel", mc.SerializeToString()])


    def plot(self,plotnames):
        plt.figure()
        plt.clf()
        legend = []
        for name in plotnames:
            legend.append(name+" fit")

            t0 = self.ts[0]
            tf = self.ts[-1]
            dt = tf-t0

            t0 -= 0.2*dt
            tf += 0.2*dt
            ts_ = np.linspace(t0,tf,500)
            plt.plot(ts_,[self.fits[name].evaluate(t) for t in ts_])

            legend.append(name)
            plt.plot(self.ts,self.fits[name].xs,'--')

        plt.title("fourier fits")
        plt.xlabel('time')
        plt.legend(legend)
        plt.grid()

    def setPhase(self,phase):
        maxPoly = max([max(self.fits[name].polyOrder) for name in self.fits])
        maxCos = max([max(self.fits[name].cosOrder) for name in self.fits])
        maxSin = max([max(self.fits[name].sinOrder) for name in self.fits])

        # all bases
        polyBases = [phase**po for po in range(maxPoly+1)]
        cosBases = [C.cos(co*phase) for co in range(maxCos+1)]
        sinBases = [C.sin(si*phase) for si in range(maxSin+1)]

        self.fitsWithPhase = {}
        for name in self.fits:
            fit = self.fits[name]

            polys = [coeff*polyBases[order] for coeff,order in zip(fit.polyCoeffs, fit.polyOrder)]
            coses = [coeff*cosBases[order]  for coeff,order in zip(fit.cosCoeffs,  fit.cosOrder)]
            sins  = [coeff*sinBases[order]  for coeff,order in zip(fit.sinCoeffs,  fit.sinOrder)]

            self.fitsWithPhase[name] = sum(polys+coses+sins)

    def lookup(self,name):
        if not hasattr(self,'fitsWithPhase'):
            raise ValueError('need to call setPhase first')
        return self.fitsWithPhase[name]

    def __getitem__(self,name):
        return self.lookup(name)

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
    filename = "data/carousel_opt"

    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    # load saved trajectory
    loadfile = filename+".dat"
    print "loading saved trajectory: "+loadfile
    f=open(loadfile,'r')
    traj = pickle.load(f)
    f.close()

    # fit everything
    xyzOrders = {'poly':[0],
                 'sin':range(1,6),
                 'cos':range(1,6)}
    xyzDotOrders = {'poly':[0],
                    'sin':range(1,8),
                    'cos':range(1,8)}
    dcmOrders = {'poly':[0],
                 'sin':range(1,10),
                 'cos':range(1,10)}
    omegaOrders = {'poly':[0],
                   'sin':range(1,7),
                   'cos':range(1,7)}
    controlSurfaceOrders = {'poly':[0],
                            'sin':range(1,5),
                            'cos':range(1,5)}
    deltaOrders = {'poly':[0,1],
                   'sin':range(1,5),
                   'cos':range(1,5)}
    ddeltaOrders = {'poly':[0],
                    'sin':range(1,5),
                    'cos':range(1,5)}
    orderMap = {'x':xyzOrders,
                'y':xyzOrders,
                'z':xyzOrders,
                'r':xyzOrders,
                'dx':xyzDotOrders,
                'dy':xyzDotOrders,
                'dz':xyzDotOrders,
                'dr':xyzDotOrders,
                'w1':omegaOrders,
                'w2':omegaOrders,
                'w3':omegaOrders,
                'e11':dcmOrders,
                'e12':dcmOrders,
                'e13':dcmOrders,
                'e21':dcmOrders,
                'e22':dcmOrders,
                'e23':dcmOrders,
                'e31':dcmOrders,
                'e32':dcmOrders,
                'e33':dcmOrders,
                'aileron':controlSurfaceOrders,
                'elevator':controlSurfaceOrders,
                'delta':deltaOrders,
                'ddelta':ddeltaOrders
                }
#    for name in traj.dvMap._xNames:
#        orderMap[name] =
    trajFit = TrajFit(orderMap,traj)

    # send kite protos
#    kiteProtos = []
#    for k in range(0,ts.size):
#        kiteProtos.append( npToKiteProto(traj.x,k,zt=0,kiteAlpha=0.2,lineAlpha=0.2) )
#        kiteProtos.append( fitToKiteProto(fits,k,zt=0, dz=1) )
#    mc = rawe.kite_pb2.MultiCarousel()
#    mc.horizon.extend(list(kiteProtos))
#    publisher.send_multipart(["multi-carousel", mc.SerializeToString()])

    savefile = filename+"_fourier.dat"
    print "saving fourier coeffs: "+savefile
    f=open(savefile,'w')
    pickle.dump(trajFit,f)
    f.close()

    def mkPlots():
        trajFit.plot(['x','y','z'])
        trajFit.plot(['r'])
        trajFit.plot(['dx','dy','dz'])
        trajFit.plot(['dr'])
        trajFit.plot(['e11','e12','e13','e21','e22','e23','e31','e32','e33'])
        trajFit.plot(['w1','w2','w3'])
        if 'delta' in traj.dvMap._xNames:
            trajFit.plot(['delta'])
        if 'ddelta' in traj.dvMap._xNames:
            trajFit.plot(['ddelta'])
        plt.show()
    mkPlots()
