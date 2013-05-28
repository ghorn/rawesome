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

from collections import Counter

def getRepeated(ns):
    c = Counter()
    for n in ns:
        c[n] += 1

    nonUnique = []
    for n,k in c.items():
        if k>1:
            nonUnique.append(n)
    return nonUnique

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
        assert len(names)==length
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
            if timestep is None:
                for timestep in range(0,self.nSteps):
                    self.dvmapSet(name,val,timestep)
                return
            # set state or action for one timestep
            val0 = self.dvmap[name][timestep]
            # warn if value being overwritten
            if (val0 is not None) and (not quiet):
                print "WARNING: "+self.descriptor+" value for \""+name+"\" at timestep "+str(timestep)+" being changed from "+str(val0)+" to "+str(val)
            self.dvmap[name][timestep] = val

        # set param
        elif name in self.pNames:
            if timestep is not None:
                raise ValueError('Can\'t set a parameter at a specific timestep')
            val0 = self.dvmap[name]
            # error if value being overwritten
            if val0 is not None:
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

    def getTimestepsFromDvs(self,dvs):
        self.lookup
        ret = {}
        ts = 0
        nx = len(self.xNames)
        nu = len(self.uNames)
        nxu = nx+nu

        xus = dvs[:self.nSteps*nxu].reshape([nxu,self.nSteps])
        p = dvs[self.nSteps*nxu:]
        
        x = []
        u = []
        for ts in range(0,self.nSteps):
            x.append(xus[:nx,ts])
            u.append(xus[nx:,ts])
        return (x,u,p)

    def devectorize(self,xup):
        ret = {}
        n = 0
        for name in self.xuNames():
            ret[name]=xup[n*self.nSteps:(n+1)*self.nSteps]
            n = n+1
        for k,name in enumerate(self.pNames):
            ret[name]=xup[n*self.nSteps+k].at(0)
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
                if val is None:
                    if name in missing:
                        missing[name].append(ts)
                    else:
                        missing[name]=[ts]
                        
        for name in self.pNames:
            if self.dvmap[name] is None:
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


    def lookup(self,name,timestep=None):
        # get state or action
        if name in self.xuNames():
            # set state or action for all timesteps
            if timestep is None:
                return [self.dvmap[name][ts] for ts in range(0,self.nSteps)]
            return self.dvmap[name][timestep]

        # get param
        elif name in self.pNames:
            if timestep is not None:
                raise ValueError('Can\'t lookup a parameter at a specific timestep')
            return self.dvmap[name]

        # error if name not in x/u/p
        else:
            raise ValueError("unrecognized variable name \""+name+"\"")
