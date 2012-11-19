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
        self.trajData.nv = ocp.getNV()
        
        self._setupOutputFuns(ocp)
        if dvs is not None:
            self.setDvs(ocp,dvs)
        
    def _setupOutputFuns(self,ocp):
        # if it's a Dae output
        self._outputFuns  = {} # no algebraic states
        self._outputFunsZ = {} # has algebraic states
        
        for name in self.trajData.outputNames:
            # try to make a function without any algebraic states
            f = C.SXFunction([ocp.dae.xVec(),ocp.dae.uVec(),ocp.dae.pVec()],
                             [ocp.dae.output(name)]
                             )
            f.init()
            if len(f.getFree()) == 0:
                self._outputFuns[name] = f
            else:
                # darn, it has algebraic states
                f = C.SXFunction([ocp.dae.xVec(),ocp.dae.zVec(),ocp.dae.uVec(),ocp.dae.pVec()],
                                 [ocp.dae.output(name)]
                                 )
                f.init()
                self._outputFunsZ[name] = f

    def setDvs(self,ocp,dvs):
        opt = ocp.devectorize(dvs)

        self.trajData.dvs = dvs
        self.trajData.collMap = opt
        self.trajData.x = {}
        self.trajData.z = {}
        self.trajData.u = {}
        self.trajData.parameters = {}
        self.trajData.outputs = {}
        self.trajData.outputsZ = {}

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

        # get outputs with no algebraic states
        for name,f in self._outputFuns.items():
            y = []
            for timestepIdx in range(opt._nk):
                for nicpIdx in range(opt._nicp):
                    for degIdx in range(opt._deg+1):
                        f.setInput(opt.xVec(timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx),0)
                        f.setInput(opt.uVec(timestepIdx),1)
                        f.setInput(opt.pVec(),2)
                        f.evaluate()
                        y.append(float(f.output(0))) # doesn't allow for vector/matrix outputs
                    y.append(numpy.nan)
            self.trajData.outputs[name] = numpy.array(y)

        # get outputs with algebraic states
        for name,f in self._outputFunsZ.items():
            y = []
            for timestepIdx in range(opt._nk):
                for nicpIdx in range(opt._nicp):
                    for degIdx in range(1,opt._deg+1):
                        f.setInput(opt.xVec(timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx),0)
                        f.setInput(opt.zVec(timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx),1)
                        f.setInput(opt.uVec(timestepIdx),2)
                        f.setInput(opt.pVec(),3)
                        f.evaluate()
                        y.append(float(f.output(0))) # doesn't allow for vector/matrix outputs
                    y.append(numpy.nan)
            self.trajData.outputsZ[name] = numpy.array(y)
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
