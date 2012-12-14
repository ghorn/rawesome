import numpy
import pickle
from collocation import Coll
import casadi as C
from trajectoryData import TrajectoryData

# non-pickleable thing which can turn design variables into x/z/u/p/outputs for TrajectoryData to use
class Trajectory(object):
    def __init__(self,ocp,dvs=None):
        assert isinstance(ocp,Coll)

        self.trajData = TrajectoryData()
        self.trajData.xNames = ocp.dae.xNames()
        self.trajData.zNames = ocp.dae.zNames()
        self.trajData.uNames = ocp.dae.uNames()
        self.trajData.pNames = ocp.dae.pNames()
        self.trajData.outputNames = ocp.dae.outputNames()
        self.trajData._outputNames0 = ocp._outputNames0
        self.trajData.nv = ocp.getNV()
        
        self._outputNamesNot0 = [ons for ons in self.trajData.outputNames
                                 if ons not in self.trajData._outputNames0]
        self._setupOutputFuns(ocp)
        
        if dvs is not None:
            self.setDvs(ocp,dvs)
        
    def _setupOutputFuns(self,ocp):
        outs = []
        for timestepIdx in range(0,ocp.nk):
            for nicpIdx in range(0,ocp.nicp):
                # deg 0:
                for name in self.trajData._outputNames0:
                    outs.append(ocp(name,timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=0))
                for name in self._outputNamesNot0:
                    # dummies to get the size right
                    outs.append(0*ocp(name,timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=1))
                # deg > 0:
                for degIdx in range(1,ocp.deg+1):
                    for name in self.trajData.outputNames:
                        outs.append(ocp(name,timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx))
        self.outputFun = C.MXFunction([ocp._dvMap.vec],outs)
        self.outputFun.init()
        
    def setDvs(self,ocp,dvs):
        opt = ocp.devectorize(dvs)

        self.trajData.dvs = dvs
        self.trajData.collMap = opt
        self.trajData.x = {}
        self.trajData.z = {}
        self.trajData.u = {}
        self.trajData.parameters = {}
        self.trajData.outputs = {}

        # make time grids
        tgrid = ocp.mkTimeGrid(dvs)
        self.trajData.tgridX = []
        self.trajData.tgridZ = []
        self.trajData.tgridU = []
        for timestepIdx in range(opt._nk):
            self.trajData.tgridU.append(tgrid[timestepIdx][0][0])
            for nicpIdx in range(opt._nicp):
                for degIdx in range(opt._deg+1):
                    self.trajData.tgridX.append(tgrid[timestepIdx][nicpIdx][degIdx])
                    if degIdx > 0:
                        self.trajData.tgridZ.append(tgrid[timestepIdx][nicpIdx][degIdx])
                # one more nan
                self.trajData.tgridX.append(numpy.nan)
                self.trajData.tgridZ.append(numpy.nan)
        self.trajData.tgridX.append(tgrid[opt._nk][0][0])
        self.trajData.tgrid = tgrid

        # get x/z/u/p
        # x
        for name in self.trajData.xNames:
            self.trajData.x[name] = []
            for timestepIdx in range(opt._nk):
                for nicpIdx in range(opt._nicp):
                    for degIdx in range(opt._deg+1):
                        self.trajData.x[name].append(opt.lookup(name,timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx))
                    self.trajData.x[name].append(numpy.nan)
            self.trajData.x[name].append(opt.lookup(name,-1,nicpIdx=0,degIdx=0))

        # z
        for name in self.trajData.zNames:
            self.trajData.z[name] = []
            for timestepIdx in range(opt._nk):
                for nicpIdx in range(opt._nicp):
                    for degIdx in range(1,opt._deg+1):
                        self.trajData.z[name].append(opt.lookup(name,timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx))
                    self.trajData.z[name].append(numpy.nan)

        # u
        for name in self.trajData.uNames:
            self.trajData.u[name] = []
            for timestepIdx in range(opt._nk):
                self.trajData.u[name].append(opt.lookup(name,timestepIdx))

        # p
        for name in self.trajData.pNames:
            self.trajData.parameters[name] = opt.lookup(name)

        # outputs
        for name in self.trajData.outputNames:
            self.trajData.outputs[name] = []
        self.outputFun.setInput(dvs,0)
        self.outputFun.evaluate()

        def grab(x):
            if x.shape == (1,1):
                return numpy.array(x[0,0])
            else:
                return x

        k = 0
        for timestepIdx in range(0,ocp.nk):
            for nicpIdx in range(0,ocp.nicp):
                # deg 0:
                for name in self.trajData._outputNames0:
                    self.trajData.outputs[name].append(grab(numpy.array(self.outputFun.output(k))))
                    k += 1
                for name in self._outputNamesNot0:
                    self.trajData.outputs[name].append(grab(numpy.nan*numpy.array(self.outputFun.output(k))))
                    k += 1
                # deg > 0:
                for degIdx in range(1,ocp.deg+1):
                    for name in self.trajData.outputNames:
                        self.trajData.outputs[name].append(grab(numpy.array(self.outputFun.output(k))))
                        k += 1
                # one more nan
                for name in self.trajData.outputNames:
                    self.trajData.outputs[name].append(numpy.nan*self.trajData.outputs[name][1])

        return self.trajData

    def plot(self,names,**kwargs):
        self.trajData.plot(names,**kwargs)
        
    def subplot(self,names,**kwargs):
        self.trajData.subplot(names,**kwargs)

    def save(self,filename):
        assert isinstance(filename,str), "filename must be a string"

        f=open(filename,'w')
        pickle.dump(self.trajData,f)
        f.close()
