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

    def toMatlab(self):
        import numpy as np

        def sanitizeName(n):
            for char in ['\'','\"',' ','(',')','-','/','\\']:
                n = n.replace(char,"_")
            if not (n[0].isalpha()):
                n = "BADNAME_"+n
            return n

        def matlabName(n):
            return n.replace("'","''")
        
        ret = ["traj = struct();"]
        
        for name in self.x:
            n = sanitizeName(name)
            ret.append("traj."+n+".name = '"+matlabName(name)+"';")
            ret.append("traj."+n+".time = "+str(self.tgridX)+";")
            ret.append("traj."+n+".data = "+str(self.x[name])+";")
    
        for name in self.u:
            n = sanitizeName(name)
            ret.append("traj."+n+".name = '"+matlabName(name)+"';")
            ret.append("traj."+n+".time = "+str(self.tgridU)+";")
            ret.append("traj."+n+".data = "+str(self.u[name])+";")

        for name in self.z:
            n = sanitizeName(name)
            ret.append("traj."+n+".name = '"+matlabName(name)+"';")
            ret.append("traj."+n+".time = "+str(self.tgridZ)+";")
            ret.append("traj."+n+".data = "+str(self.z[name])+";")

        for name in self.outputs:
            out = self.outputs[name]
            for k in range(out.size):
                try:
                    out[k] = out[k][0][0]
                except TypeError:
                    pass
            n = sanitizeName(name)
            ret.append("traj."+n+".name = '"+matlabName(name)+"';")
            ret.append("traj."+n+".time = "+str(self.tgridX[:-1])+";")
#            ret.append("traj."+n+".data = "+str(self.outputs[name].flatten())+";")
            ret.append("traj."+n+".data = "+str(out.tolist())+";")

        return "\n".join(ret)
                       
    def saveAsMatlab(self,filename):
        out = self.toMatlab()
        f=open(filename,'w')
        f.write(out)
        f.close()
        
    
if __name__=='__main__':
    import pickle
    
    filename = "data/crosswind_opt.dat"
    print "loading "+filename
    f=open(filename,'r')
    traj = pickle.load(f)
    f.close()

    traj.saveAsMatlab('data/crosswind_opt.m')
