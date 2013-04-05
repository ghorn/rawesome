import casadi as C

import exportOcp

class Ocp(object):
    def __init__(self, dae, N=None, ts=None):
        dae.assertNoFreeParams()
        self._dae = dae
        if N is None or ts is None:
            raise Exception('please initialize Ocp with Ocp(dae, N=.., ts=..)')
        assert type(N) is int, "N must be an int, got type: "+str(type(N))
        assert type(ts) in [int,float], "ts must be an int or float, got type: "+str(type(ts))
        self._nk = N
        self._ts = float(ts)
#        self._bndmap = {}
#        self._bndmapStart = {}
#        self._bndmapEnd = {}

        self._constraints = []
        self._constraintsStart = []
        self._constraintsEnd = []

    # stuff inherited from dae
    def __getitem__(self,name):
        return self._dae[name]
    def __contains__(self,name):
        if not isinstance(name,str):
            raise KeyError('key must be a string')
        return name in self._dae
    def ddt(self,name):
        return self._dae.ddt(name)


#    def bound(self, name, (lb,ub), when=None):
#        assert name in dae, "unrecognized name \""+name+"\""
#        def blah(bmap,msg):
#            if name in bmap:
#                msg = '"'+name+'" has already been bound'+msg+', old bounds: '+\
#                    str(bmap[name])+', new bounds: '+str((lb,ub))
#                raise Exception(msg)
#            else:
#                bmap[name] = (lb,ub)
#
#        if when is None:
#            blah(self._bndmap, '')
#        elif when is "AT_END":
#            blah(self._bndmapEnd, ' AT_END')
#        elif when is "AT_START":
#            blah(self._bndmapStart, ' AT_START')
#        else:
#            raise Exception("unrecognized \"when\": "+str(when))

    def constrain(self,*args, **kwargs):
        usage = "do Ocp.constrain(x, CMP, y) or Ocp.constrain(x, CMP1, y, CMP2, z) where CMP,CMP1,CMP2 can be any of '<=', '==', or '>=' written as strings"
        if len(args) is 3:
            assert args[1] in ['>=','==','<='], usage
            self._constrainOne(*args, **kwargs)
        elif len(args) is 5:
            assert args[1] in ['>=','==','<='], usage
            assert args[3] in ['>=','==','<='], usage
            self._constrainOne(args[0], args[1], args[2], **kwargs)
            self._constrainOne(args[2], args[3], args[4], **kwargs)
        else:
            raise Exception(usage)

    def _constrainOne(self, lhs, comparison, rhs, when=None):
        if type(lhs) == C.SXMatrix:
            assert lhs.shape == (1,1), "lhs must be scalar, got matrix with shape: "+str(lhs.shape)
        else:
            assert type(lhs) in [int,float], "lhs type unrecognized: "+str(type(lhs))
        if type(rhs) == C.SXMatrix:
            assert rhs.shape == (1,1), "rhs must be scalar, got matrix with shape: "+str(rhs.shape)
        else:
            assert type(rhs) in [int,float], "rhs type unrecognized: "+str(type(rhs))

        assert comparison in ['==','<=','>='], 'comparison "'+str(comparison)+\
            '" is not "==", "<=", or ">="'
        if when is None:
            self._constraints.append((lhs,comparison,rhs))
        elif when is "AT_END":
            self._constraintsEnd.append((lhs,comparison,rhs))
        elif when is "AT_START":
            self._constraintsStart.append((lhs,comparison,rhs))
        else:
            raise Exception("\"when\" must be 'AT_START' or 'AT_END', leaving it blank means always, you put: "+str(when))

    def minimizeLsq(self, obj):
        C.makeDense(obj)
        shape = obj.shape
        assert shape[0] == 1 or shape[1] == 1, 'objective cannot be matrix, got shape: '+str(shape)
        assert not hasattr(self, '_minLsq'), 'you can only call minimizeLsq once'
        self._minLsq = obj

    def minimizeLsqEndTerm(self, obj):
        C.makeDense(obj)
        shape = obj.shape
        assert shape[0] == 1 or shape[1] == 1, 'objective cannot be matrix, got shape: '+str(shape)
        assert not hasattr(self, '_minLsqEndTerm'), 'you can only call minimizeLsqEndTerm once'
        self._minLsqEndTerm = obj

    def exportCode(self, options, qpSolver='QP_OASES'):
        print "put default options in here"
        return exportOcp.exportOcp(self, options, qpSolver)
