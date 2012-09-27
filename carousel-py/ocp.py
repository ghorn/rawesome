import numpy
from collections import Counter
import numbers

import casadi as C

def getRepeated(ns):
    c = Counter()
    for n in ns:
        c[n] += 1

    nonUnique = []
    for n,k in c.items():
        if k>1:
            nonUnique.append(n)
    return nonUnique


class Constraints():
    def __init__(self):
        self._g = []
        self._glb = []
        self._gub = []
        
    def add(self,lhs,comparison,rhs):
        #print "\n\nadding constraint\nlhs: "+str(lhs)+"\ncomparison: "+comparison+"\nrhs: "+str(rhs)
        if comparison=="==":
            g = lhs - rhs
            self._g.append(g)
            self._glb.append(numpy.zeros(g.size()))
            self._gub.append(numpy.zeros(g.size()))
        elif comparison=="<=":
            g = lhs - rhs
            self._g.append(g)
            self._glb.append(-numpy.inf*numpy.ones(g.size()))
            self._gub.append(numpy.zeros(g.size()))
        elif comparison==">=":
            g = rhs - lhs
            self._g.append(g)
            self._glb.append(-numpy.inf*numpy.ones(g.size()))
            self._gub.append(numpy.zeros(g.size()))
        else:
            raise ValueError('Did not recognize comparison \"'+str(comparison)+'\"')

    def addDynamicsConstraints(self,integrator,states,actions,params=None):
        nSteps = states.size2()
        if nSteps != actions.size2():
            raise ValueError("actions and states have different number of steps")

        for k in range(0,nSteps-1):
            u = actions[:,k]
            if params != None: # params are appended to control inputs
                u = C.veccat([u,params])
            xk   = states[:,k]
            xkp1 = states[:,k+1]
            self.add(integrator.call([xk,u])[C.INTEGRATOR_XF],'==',xkp1)

    def getG(self):
        return C.veccat(self._g)
    def getLb(self):
        return C.veccat(self._glb)
    def getUb(self):
        return C.veccat(self._gub)

class DesignVarMap():
    descriptor = ""
    def __init__(self, xNames, uNames, pNames, nSteps):
        r = getRepeated(xNames+uNames+pNames)
        if len(r)>0:
            raise ValueError("there are redundant names in the OCP: "+str(r))

        self.nSteps = nSteps
        self.xNames = xNames
        self.uNames = uNames
        self.pNames = pNames
        
        self.dvmap = {}
        
        for name in self.xuNames():
            self.dvmap[name] = [None for k in range(0,self.nSteps)]
        for name in self.pNames:
            self.dvmap[name] = None

    def xuNames(self):
        return self.xNames+self.uNames

    def _dvmapSetVec(self,val,names,**kwargs):
        if isinstance(val,list):
            length = len(val)
        elif hasattr(val,'size'):
            length = val.size()
        else:
            raise ValueError("can't figure out how long "+str(val)+" is")
        assert(len(names)==length)
        for k,name in enumerate(names):
            self.dvmapSet(name,val[k],**kwargs)
        
    def setXVec(self,val,**kwargs):
        self._dvmapSetVec(val,self.xNames,**kwargs)

    def setUVec(self,val,**kwargs):
        self._dvmapSetVec(val,self.uNames,**kwargs)

    def dvmapSet(self,name,val,timestep=None,quiet=False):
        # set state or action
        if name in self.xuNames():
            # set state or action for all timesteps
            if timestep==None:
                for timestep in range(0,self.nSteps):
                    self.dvmapSet(name,val,timestep)
                return
            # set state or action for one timestep
            val0 = self.dvmap[name][timestep]
            # warn if value being overwritten
            if val0 != None and not quiet:
                print "WARNING: "+self.descriptor+" value for \""+name+"\" at timestep "+str(timestep)+" being changed from "+str(val0)+" to "+str(val)
            self.dvmap[name][timestep] = val

        # set param
        elif name in self.pNames:
            if timestep!=None:
                raise ValueError('Can\'t set a parameter at a specific timestep')
            val0 = self.dvmap[name]
            # error if value being overwritten
            if val0 != None:
                raise ValueError(self.descriptor+" value for parameter \""+name+"\" being changed from "+str(val0)+" to "+str(val))
            self.dvmap[name] = val

        # error if name not in x/u/p
        else:
            raise ValueError("unrecognized variable name \""+name+"\"")

    def vectorize(self):
        # make sure all bounds are set
        self._assertAllValuesSet()
        
        # concatenate then unzip bounds
        return self._concatValues()

    def devectorize(self,xup):
        ret = {}
        n = 0
        for name in self.xuNames():
            ret[name]=xup[n*self.nSteps:(n+1)*self.nSteps]
            n = n+1
        for k,name in enumerate(self.pNames):
            ret[name]=xup[n*self.nSteps+k]
        return ret

    def _concatValues(self):
        xuVals = [self.dvmap[name] for name in self.xuNames()]
         
        import itertools
        chain = itertools.chain(*xuVals)
        return list(chain)+[self.dvmap[name] for name in self.pNames]
    
    def _assertAllValuesSet(self): # make sure all bounds are set
        # populate dictionary of missing values
        missing = {}
        for name in self.xuNames():
            for ts,val in enumerate(self.dvmap[name]):
                if val==None:
                    if name in missing:
                        missing[name].append(ts)
                    else:
                        missing[name]=[ts]
                        
        for name in self.pNames:
            if self.dvmap[name]==None:
                missing[name] = True

        # if dictionary is not empty, raise error
        errs = []
        for name in self.xuNames()+self.pNames:
            if name in missing:
                if isinstance(missing[name],list):
                    errs.append(self.descriptor+" missing for: state/action \""+name+"\" at timesteps "+str(missing[name]))
                else:
                    errs.append(self.descriptor+" missing for: parameter \""+name+"\"")
        if len(errs)>0:
            raise ValueError("\n"+"\n".join(errs))
                

class Bounds(DesignVarMap):
    descriptor = "bound"
        
    def setBound(self,name,val,**kwargs):
        assert(isinstance(name,str))
        assert(isinstance(val,tuple))
        assert(len(val)==2)
        assert(isinstance(val[0],numbers.Real))
        assert(isinstance(val[1],numbers.Real))
        self.dvmapSet(name,val,**kwargs)

    def get(self):
        return zip(*DesignVarMap.vectorize(self))

class InitialGuess(DesignVarMap):
    descriptor = "initial guess"

    def setGuess(self,name,val,**kwargs):
        assert(isinstance(name,str))
        assert(isinstance(val,numbers.Real))
        self.dvmapSet(name,val,**kwargs)
