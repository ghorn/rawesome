import numpy
from collections import Counter

import casadi as C


class Constraints():
    def __init__(self):
        self._g = []
        self._glb = []
        self._gub = []
        
    def add(self,lhs,comparison,rhs):
        if comparison=="==":
            g = lhs - rhs
            self._g.append(g)
            self._glb.append(numpy.zeros(g.size()))
            self._gub.append(numpy.zeros(g.size()))
        elif comparison=="<=":
            g = lhs - rhs
            self._g.append(g)
            self._glb.append(-numpy.inf*numpy.ones(h.size()))
            self._gub.append(numpy.zeros(g.size()))
        elif comparison==">=":                                         
            g = rhs - lhs
            self._g.append(g)
            self._glb.append(-numpy.inf*numpy.ones(h.size()))
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
            integrator.call([xk,u])
            self.add(integrator.call([xk,u])[C.INTEGRATOR_XF],'==',xkp1)

    def getG(self):
        return C.veccat(self._g)
    def getLb(self):
        return C.veccat(self._glb)
    def getUb(self):
        return C.veccat(self._gub)

class Bounds():
    def __init__(self, nSteps, xNames, uNames, pNames):
        def getRepeated(ns):
            c = Counter()
            for n in ns:
                c[n] += 1

            nonUnique = []
            for n,k in c.items():
                if k>1:
                    nonUnique.append(n)
            return nonUnique

        r = getRepeated(xNames+uNames+pNames)
        if len(r)>0:
            raise ValueError("there are redundant names in the OCP: "+str(r))

        self.nSteps = nSteps
        self.xuNames = xNames+uNames
        self.pNames = pNames
        
        self.bounds = {}
        
        for name in self.xuNames:
            self.bounds[name] = [None for k in range(0,self.nSteps)]
        for name in self.pNames:
            self.bounds[name] = None

    def bound(self,name,lbub,timestep=None):
        # set state or action
        if name in self.xuNames:
            # set state or action for all timesteps
            if timestep==None:
                for timestep in range(0,self.nSteps):
                    self.bound(name,lbub,timestep)
                return
            # set state or action for one timestep
            bnd0 = self.bounds[name][timestep]
            # warn if bound being overwritten
            if bnd0 != None:
                print "WARNING: bound for \""+name+"\" at timestep "+str(timestep)+" being changed from "+str(bnd0)+" to "+str(lbub)
            self.bounds[name][timestep] = lbub

        # set param
        elif name in self.pNames:
            if timestep!=None:
                raise ValueError('Can\'t bound a parameter at a specific timestep')
            bnd0 = self.bounds[name]
            # error if bound being overwritten
            if bnd0 != None:
                raise ValueError("bound for parameter \""+name+"\" being changed from "+str(bnd0)+" to "+str(lbub))
            self.bounds[name] = lbub

        # error if name not in x/u/p
        else:
            raise ValueError("unrecognized variable name \""+name+"\"")

    def get(self):
        # make sure all bounds are set
        self._checkBnds()
        
        # concatenate then unzip bounds
        return zip(*self._concatBnds())

    def _concatBnds(self):
        xuBounds = [self.bounds[name] for name in self.xuNames]
         
        import itertools
        chain = itertools.chain(*xuBounds)
        return list(chain)+[self.bounds[name] for name in self.pNames]
    
    def _checkBnds(self): # make sure all bounds are set
        # populate dictionary of missing values
        missing = {}
        for name in self.xuNames:
            for ts,val in enumerate(self.bounds[name]):
                if val==None:
                    if name in missing:
                        missing[name].append(ts)
                    else:
                        missing[name]=[ts]
                        
        for name in self.pNames:
            if self.bounds[name]==None:
                missing[name] = True

        # if dictionary is not empty, raise error
        errs = []
        for name in self.xuNames+self.pNames:
            if name in missing:
                if isinstance(missing[name],list):
                    errs.append("missing bound for: state/action \""+name+"\" at timesteps "+str(missing[name]))
                else:
                    errs.append("missing bound for: parameter \""+name+"\"")
        if len(errs)>0:
            raise ValueError("\n"+"\n".join(errs))
                

                
