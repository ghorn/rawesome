import matplotlib.pyplot as plt
import numpy
import pickle

# pickleable thing which stores x/z/u/p/outputs and can plot them
class TrajectoryData(object):
    xNames = None
    zNames = None
    uNames = None
    pNames = None
    outputNames = None
    nv = None
    
    tgrid = None
    tgridZOutput = None
    xzu = None
    parameters = None
    outputs = None
    outputsZ = None

    def subplot(self,names,title=None):
        assert isinstance(names,list)
        
        plt.figure()
        plt.clf()
        n = len(names)
        for k,name in enumerate(names):
            plt.subplot(n,1,k+1)
            if k==0:
                self._plot(name,title)
            else:
                self._plot(name,None)

    def plot(self,names,title=None):
        plt.figure()
        plt.clf()
        self._plot(names,title)

    def _plot(self,names,title):
        if isinstance(names,str):
            names = [names]
        assert isinstance(names,list)
        legend = []
        for name in names:
            assert isinstance(name,str)
            legend.append(name)

            # if it's a simple state or action
            if name in self.xzu:
                plt.plot(self.tgrid,self.xzu[name])

            # if it's a dae output with NO algebraic states
            elif name in self.outputs:
                plt.plot(self.tgrid,self.outputs[name])

            # if it's a dae output WITH algebraic states
            elif name in self.outputsZ:
                plt.plot(self.tgridZOutput,self.outputsZ[name])

            # throw error on parameter
            elif name in self.parameters:
                raise ValueError("can't plot a parameter (\""+name+"\")")

            # throw error on unrecognized
            else:
                raise ValueError("unrecognized name \""+name+"\"")
        
        if title is not None:
            assert isinstance(title,str), "title must be a string"
            plt.title(title)
        plt.xlabel('time')
        plt.legend(legend)
        plt.grid()


# non-pickleable thing which can turn design variables into x/z/u/p/outputs for TrajectoryData to use
class Trajectory(object):
    def __init__(self,ocp):
        from collocation import Coll
        assert isinstance(ocp,Coll)

        self.trajData = TrajectoryData()
        self.trajData.xNames = ocp.dae.xNames()
        self.trajData.zNames = ocp.dae.zNames()
        self.trajData.uNames = ocp.dae.uNames()
        self.trajData.pNames = ocp.dae.pNames()
        self.trajData.outputNames = ocp.dae.outputNames()
        self.trajData.nv = ocp.getNV()
        
        self._setupOutputFuns(ocp)
        
    def _setupOutputFuns(self,ocp):
        import casadi as C

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

        self.trajData.xzu = {}
        self.trajData.parameters = {}
        self.trajData.outputs = {}
        self.trajData.outputsZ = {}
        
        self.trajData.tgrid = opt['tgrid']
        self.trajData.tgridZOutput = opt['tgrid'][:-1]

        # get x/u/z
        for name in self.trajData.xNames+self.trajData.uNames+self.trajData.zNames:
            self.trajData.xzu[name] = opt['vardict'][name]

        # get parameters
        for name in self.trajData.pNames:
            self.trajData.parameters[name] = opt['vardict'][name]

        # get outputs with no algebraic states
        for name,f in self._outputFuns.items():
            y = []
            for k,t in enumerate(self.trajData.tgrid):
                f.setInput(opt['x'][:,k],0)
                f.setInput(opt['u'][:,k],1)
                f.setInput(opt['p'],2)
                f.evaluate()
                y.append(float(f.output(0))) # doesn't allow for vector/matrix outputs
            self.trajData.outputs[name] = numpy.array(y)

        # get outputs with algebraic states
        for name,f in self._outputFunsZ.items():
            y = []
            for k,t in enumerate(self.trajData.tgridZOutput):
                f.setInput(opt['x'][:,k],0)
                f.setInput(opt['zPlt'][:,k],1)
                f.setInput(opt['u'][:,k],2)
                f.setInput(opt['p'],3)
                f.evaluate()
                y.append(float(f.output(0))) # doesn't allow for vector/matrix outputs
            self.trajData.outputsZ[name] = numpy.array(y)
        return self.trajData

    def plot(self,*args,**kwargs):
        self.trajData.plot(*args,**kwargs)
        
    def subplot(self,*args,**kwargs):
        self.trajData.subplot(*args,**kwargs)

    def save(self,filename):
        assert isinstance(filename,str), "filename must be a string"

        f=open(filename,'w')
        pickle.dump(opt,f)
        f.close()
