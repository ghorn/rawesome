import numpy as np

class CollMap(object):
    def __init__(self,ocp,name):
        # grab everything needed from the ocp
        self._xNames = ocp.dae.xNames()
        self._zNames = ocp.dae.zNames()
        self._uNames = ocp.dae.uNames()
        self._pNames = ocp.dae.pNames()
        self._nk = ocp.nk
        self._nicp = ocp.nicp
        self._deg = ocp.deg
        self._collPoly = ocp.collPoly

        assert isinstance(name,str)
        self._name = name
        self._initMap()

    def _initMap(self):
        self._xMap = {}
        self._zMap = {}
        self._uMap = {}
        self._pMap = {}
        for name in self._xNames:
            self._xMap[name] = np.resize(np.array([None]),(self._nk+1,self._nicp,self._deg+1))
        for name in self._zNames:
            self._zMap[name] = np.resize(np.array([None]),(self._nk,self._nicp,self._deg+1))
        for name in self._uNames:
            self._uMap[name] = np.resize(np.array([None]),(self._nk))
        for name in self._pNames:
            self._pMap[name] = None

    def lookup(self,name,timestep=None,nicpIdx=None,degIdx=None):
        return self._lookupOrSet(name,timestep,nicpIdx,degIdx)
        
    def setVal(self,name,val,timestep=None,nicpIdx=None,degIdx=None,quiet=False,force=False):
        self._lookupOrSet(name,timestep,nicpIdx,degIdx,setVal=val,quiet=quiet,force=force)
        
    def _lookupOrSet(self,name,timestep,nicpIdx,degIdx,setVal=None,quiet=False,force=False):
        assert isinstance(name,str), "lookup key must be a string in "+self._name+" map"
        
        if name in self._xMap:
            assert timestep is not None, "must give timestep for differential state lookup ("+self._name+")"
            if nicpIdx is None:
                nicpIdx = 0
            if degIdx is None:
                degIdx = 0
            assert timestep < (self._nk+1), \
                   "timestep: "+str(timestep)+" out of range in "+self._name+" map (nk: "+str(nk)+")"
            assert degIdx >=0 and degIdx < (self._deg+1), \
                   "degIdx: "+str(deg)+" out of range in "+self._name+" map (deg: "+str(self._deg)+")"
            if timestep is self._nk:
                assert nicpIdx==0 and degIdx==0,"last timestep is only defined at nicpIdx=0,degIdx=0"
            if setVal is None:
                return self._xMap[name][timestep][nicpIdx][degIdx]
            else:
                oldval = self._xMap[name][timestep][nicpIdx][degIdx]
                if (not quiet) and (oldval is not None):
                    print "WARNING: changing \"%s\" %s at timestep %d from %s to %s" % (name,self._name,timestep,str(oldval),str(setVal))
                self._xMap[name][timestep][nicpIdx][degIdx] = setVal

        elif name in self._zMap:
            assert timestep is not None, "must give timestep for algebraic state lookup ("+self._name+")"
            if nicpIdx is None:
                nicpIdx = 0
            assert degIdx is not None, "must set degIdx for algebraic state map ("+self._name+")"
            assert timestep < self._nk, \
                   "timestep: "+str(timestep)+" out of range in "+self._name+" map (nk: "+str(nk)+")"
            assert degIdx >=0 and degIdx < (self._deg+1), \
                   "degIdx: "+str(degIdx)+" out of range in "+self._name+" map (deg: "+str(self._deg)+")"
            if setVal is None:
                return self._zMap[name][timestep][nicpIdx][degIdx]
            else:
                oldval = self._zMap[name][timestep][nicpIdx][degIdx]
                if (quiet is False) and (oldval is not None):
                    print "WARNING: changing \"%s\" %s at timestep %d from %s to %s" % (name,self._name,timestep,str(oldval),str(setVal))
                self._zMap[name][timestep][nicpIdx][degIdx] = setVal

        elif name in self._uMap:
            assert timestep is not None, "must give timestep for control input lookup ("+self._name+")"
            assert nicpIdx is None, "nicpIdx invalid for control input ("+self._name+")"
            assert degIdx is None, "degIdx invalid for control input ("+self._name+")"
            assert timestep < self._nk, \
                   "timestep: "+str(timestep)+" out of range in "+self._name+" map (nk: "+str(nk)+")"
            if setVal is None:
                return self._uMap[name][timestep]
            else:
                oldval = self._uMap[name][timestep]
                if (quiet is False) and (oldval is not None):
                    print "WARNING: changing \"%s\" %s at timestep %d from %s to %s" % (name,self._name,timestep,str(oldval),str(setVal))
                self._uMap[name][timestep] = setVal

        elif name in self._pMap:
            assert timestep is None, "timestep invalid for parameter lookup ("+self._name+")"
            assert nicpIdx is None, "nicpIdx invalid for parameter lookup ("+self._name+")"
            assert degIdx is None, "degIdx invalid for parameter lookup ("+self._name+")"
            if setVal is None:
                return self._pMap[name]
            else:
                oldval = self._pMap[name]
                if (force is False) and (oldval is not None):
                    raise ValueError("can't change \""+name+"\" "+self._name+" once it's set unless you use force=True (tried to change "+str(oldval)+" to "+str(setVal))
                self._pMap[name] = setVal

        else:
            raise KeyError("couldn't find \""+name+"\" in "+self._name+" map")
