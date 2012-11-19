import matplotlib.pyplot as plt

# pickleable thing which stores x/z/u/p/outputs and can plot them
class TrajectoryData(object):
    def subplot(self,names,title=None):
        assert isinstance(names,list)
        
        fig = plt.figure()
        if title is None:
            if isinstance(names,str):
                title = names
            else:
                assert isinstance(names,list)
                if len(names) == 1:
                    title = names[0]
                else:
                    title = str(names)
        fig.canvas.set_window_title(str(title))
                    
        plt.clf()
        n = len(names)
        for k,name in enumerate(names):
            plt.subplot(n,1,k+1)
            if k==0:
                self._plot(name,title)
            else:
                self._plot(name,None)

    def plot(self,names,title=None):
        fig = plt.figure()
        if title is None:
            if isinstance(names,str):
                title = names
            else:
                assert isinstance(names,list)
                if len(names) == 1:
                    title = names[0]
                else:
                    title = str(names)
        
        fig.canvas.set_window_title(str(title))
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
            if name in self.x:
                plt.plot(self.tgridX,self.x[name])
            elif name in self.z:
                plt.plot(self.tgridZ,self.z[name])
            elif name in self.u:
                plt.plot(self.tgridU,self.u[name],drawstyle='steps-pre')

            # if it's a dae output with NO algebraic states
            elif name in self.outputs:
                plt.plot(self.tgridX[:-1],self.outputs[name])

            # if it's a dae output WITH algebraic states
            elif name in self.outputsZ:
                plt.plot(self.tgridZ,self.outputsZ[name])

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
