import matplotlib.pyplot as plt

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
