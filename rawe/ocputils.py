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

    def add(self,lhs,comparison,rhs,tag=('unnamed_constraint',None)):
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

    def addBnds(self,g,(glb,gub),tag=('unnamed_constraint',None)):
        (tagName, tagIdx) = tag

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

    def getViolations(self,g,lbg,ubg,reportThreshold=0,reportEqViolations=False):
        """
        Tests if g >= ubg + reportThreshold
                 g <= lbg - reportThreshold
        Positive reportThreshold supresses barely active bounds
        Negative reportThreshold reports not-quite-active bounds
        """
        violations = {}

        ubviols = g - ubg
        lbviols = lbg - g
        ubviolsIdx = np.where(C.logic_and(ubviols >= reportThreshold, ubg > lbg))[0]
        lbviolsIdx = np.where(C.logic_and(lbviols >= reportThreshold, ubg > lbg))[0]

        eqviolsIdx = []
        if reportEqViolations:
            eqviolsIdx = np.where(C.logic_and(C.fabs(ubviols) >= reportThreshold, ubg == lbg))[0]

        violations = {}
        for k in ubviolsIdx:
            (name,time,idx) = self._tags[k]
            viol = ('ub',(time,idx),float(ubviols[k]))#,g[k],ubg[k])
            if name not in violations:
                violations[name] = [viol]
            else:
                violations[name].append(viol)
        for k in lbviolsIdx:
            (name,time,idx) = self._tags[k]
            viol = ('lb',(time,idx),float(lbviols[k]))#,g[k],lbg[k])
            if name not in violations:
                violations[name] = [viol]
            else:
                violations[name].append(viol)
        for k in eqviolsIdx:
            (name,time,idx) = self._tags[k]
            viol = ('eq',(time,idx),float(ubviols[k]))#,g[k],lbg[k])
            if name not in violations:
                violations[name] = [viol]
            else:
                violations[name].append(viol)
        return violations

    def getViolationsStr(self,*args,**kwargs):
        viols = self.getViolations(*args,**kwargs)
        ret = []
#        from termcolor import colored
        for name in viols:
            vstrs = str(sorted(viols[name], key=lambda x: -x[2]))
#            vstrs = ["("+str(degidx)+","+str(k)+","+colored(str(val))+")" for (degidx,k,val) in sorted(viols[name], key=lambda x: -x[2])]
            blah = "constraint violation! \""+name+": "+vstrs
            if len(blah) > 150:
                blah = blah[:150] + "..."
            ret.append(blah)
        return '\n'.join(ret)

    def printViolations(self,*args,**kwargs):
        viols = self.getViolationsStr(*args,**kwargs)
        if viols != '':
            print viols

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
