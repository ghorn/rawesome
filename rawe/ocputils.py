import numpy as np
import numbers

#from models import Dae
from dvmap import DesignVarMap

import casadi as C

def setFXOptions(fun, options):
    assert isinstance(options,list)
    for intOpt in options:
        assert isinstance(intOpt,tuple)
        assert len(intOpt)==2
        assert isinstance(intOpt[0],str)
        optName,optVal = intOpt
        fun.setOption(optName, optVal)


class Constraints():
    def __init__(self):
        self._tags = []
        self._g = []
        self._glb = []
        self._gub = []
        
    def add(self,lhs,comparison,rhs,tag):
        #print "\n\nadding constraint\nlhs: "+str(lhs)+"\ncomparison: "+comparison+"\nrhs: "+str(rhs)
        if comparison=="==":
            g = lhs - rhs
            glb = np.zeros(g.size())
            gub = np.zeros(g.size())
            self.addBnds(g,(glb,gub),tag)
        elif comparison=="<=":
            g = lhs - rhs
            glb = -np.inf*np.ones(g.size())
            gub = np.zeros(g.size())
            self.addBnds(g,(glb,gub),tag)
        elif comparison==">=":
            g = rhs - lhs
            glb = -np.inf*np.ones(g.size())
            gub = np.zeros(g.size())
            self.addBnds(g,(glb,gub),tag)
        else:
            raise ValueError('Did not recognize comparison \"'+str(comparison)+'\"')

    def addBnds(self,g,(glb,gub),(tagName,tagIdx)):
        if (isinstance(glb,numbers.Real) and isinstance(gub,numbers.Real)):
            glb = np.array(glb)
            gub = np.array(gub)

        assert isinstance(glb,np.ndarray)
        assert isinstance(gub,np.ndarray)
        assert isinstance(g,C.SXMatrix) or isinstance(g,C.MX)
        assert g.size()==glb.size and g.size()==gub.size
        self._g.append(g)
        self._glb.append(glb)
        self._gub.append(gub)
        for k in range(g.size()):
            self._tags.append( (tagName,tagIdx,k) )

    def getG(self):
        return C.veccat(self._g)
    def getLb(self):
        return C.veccat(self._glb)
    def getUb(self):
        return C.veccat(self._gub)

    def getViolations(self,g,lbg,ubg,reportThreshold=0):
        """
        Tests if g >= ubg + reportThreshold
                 g <= lbg - reportThreshold
        Positive reportThreshold supresses barely active bounds
        Negative reportThreshold reports not-quite-active bounds
        """
        violations = {}

        ubviols = g - ubg
        lbviols = lbg - g
        ubviolsIdx = np.where(ubviols >= reportThreshold)[0]
        lbviolsIdx = np.where(lbviols >= reportThreshold)[0]
        violations = {}
        for k in ubviolsIdx:
            (name,time,idx) = self._tags[k]
            viol = ('ub',(time,idx),float(ubviols[k]))
            if name not in violations:
                violations[name] = [viol]
            else:
                violations[name].append(viol)
        for k in lbviolsIdx:
            (name,time,idx) = self._tags[k]
            viol = ('lb',(time,idx),float(lbviols[k]))
            if name not in violations:
                violations[name] = [viol]
            else:
                violations[name].append(viol)
        return violations

    def printViolations(self,*args,**kwargs):
#        from termcolor import colored
        viols = self.getViolations(*args,**kwargs)
        for name in viols:
            vstrs = str(sorted(viols[name], key=lambda x: -x[2]))
#            vstrs = ["("+str(degidx)+","+str(k)+","+colored(str(val))+")" for (degidx,k,val) in sorted(viols[name], key=lambda x: -x[2])]
            print "constraint violation! \""+name+": "+vstrs[:150]

class Bounds(DesignVarMap):
    descriptor = "bound"
        
    def setBound(self,name,val,**kwargs):
        assert isinstance(name,str)
        assert isinstance(val,tuple)
        assert len(val)==2
        assert isinstance(val[0],numbers.Real)
        assert isinstance(val[1],numbers.Real)
        self.dvmapSet(name,val,**kwargs)

    def get(self):
        return zip(*DesignVarMap.vectorize(self))

class InitialGuess(DesignVarMap):
    descriptor = "initial guess"

    def setGuess(self,name,val,**kwargs):
        assert isinstance(name,str)
        assert isinstance(val,numbers.Real)
        self.dvmapSet(name,val,**kwargs)

class DesignVars(DesignVarMap):
    descriptor = "design variables"
    def __init__(self, (xNames,states), (uNames,actions), (pNames,params), nSteps):
        DesignVarMap.__init__(self,xNames,uNames,pNames,nSteps)
        for k in range(0,self.nSteps):
            self.setXVec(states[:,k],timestep=k)
            self.setUVec(actions[:,k],timestep=k)
        for k,pname in enumerate(pNames):
            self.dvmapSet(pname,params[k])
