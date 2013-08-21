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

import casadi as C

import exportOcp
from ..rtIntegrator import RtIntegratorOptions
from ..utils.options import Options, OptStr, OptInt, OptBool, OptDouble

class OcpExportOptions(Options):
    def __init__(self):
        Options.__init__(self, 'OCP')
        self.add(OptStr('SPARSE_QP_SOLUTION',
                        ['CONDENSING','FULL_CONDENSING','FULL_CONDENSING_N2','SPARSE_SOLVER']))
        self.add(OptStr('QP_SOLVER',['QP_QPOASES','QP_QPDUNES','QP_FORCES']))
        self.add(OptStr('HESSIAN_APPROXIMATION',['GAUSS_NEWTON'],default='GAUSS_NEWTON'))
        # TODO Hide this and set it as default, MS
        self.add(OptStr('DISCRETIZATION_TYPE',['MULTIPLE_SHOOTING'],default='MULTIPLE_SHOOTING'))
        self.add(OptBool('GENERATE_TEST_FILE',default=False))
        self.add(OptBool('GENERATE_MAKE_FILE',default=False))
        self.add(OptBool('GENERATE_MATLAB_INTERFACE',default=False))
        self.add(OptBool('HOTSTART_QP',default=False))
        self.add(OptBool('FIX_INITIAL_STATE',default=True))
        self.add(OptDouble('LEVENBERG_MARQUARDT',default=0.0))
        self.add(OptBool('CG_USE_VARIABLE_WEIGHTING_MATRIX',default=False))
#        self.add(OptBool('CG_USE_C99',default=True))

class Ocp(object):
    def __init__(self, dae, N=None, ts=None):
        dae.assertNoFreeParams()
        self.hashPrefix = 'ocp'
        self._dae = dae
        if N is None or ts is None:
            raise Exception('please initialize Ocp with Ocp(dae, N=.., ts=..)')
        assert type(N) is int, "N must be an int, got type: "+str(type(N))
        assert type(ts) in [int,float], "ts must be an int or float, got type: "+str(type(ts))
        self._N = N
        self._ts = float(ts)

        self._ebndmap = {}
        self._ebndmapStart = {}
        self._ebndmapEnd = {}

        self._lbndmap = {}
        self._lbndmapStart = {}
        self._lbndmapEnd = {}

        self._ubndmap = {}
        self._ubndmapStart = {}
        self._ubndmapEnd = {}

        self._constraints = []
        self._constraintsStart = []
        self._constraintsEnd = []

        self._dbgMessages = []

    @property
    def N(self):
        return self._N
    @property
    def ts(self):
        return self._ts
    @property
    def dae(self):
        return self._dae

    # stuff inherited from dae
    def __getitem__(self,name):
        return self.dae[name]
    def __contains__(self,name):
        if not isinstance(name,str):
            raise KeyError('key must be a string')
        return name in self.dae
    def ddt(self,name):
        return self.dae.ddt(name)

    # debugging message
    def debug(self,msg):
        self._dbgMessages.append(msg)

    # add a linear constraint
    def _bound(self, name, bnd, upperLowerEq, when=None):
        assert name in self.dae.xNames()+self.dae.uNames(), "you can't bound "+name+\
            " because only differential state and control bounds are supported"

        if upperLowerEq == 'lower':
            maps = {'ALWAYS':self._lbndmap,
                    'AT_START':self._lbndmapStart,
                    'AT_END':self._lbndmapEnd,
                    }
        elif upperLowerEq == 'upper':
            maps = {'ALWAYS':self._ubndmap,
                    'AT_START':self._ubndmapStart,
                    'AT_END':self._ubndmapEnd,
                    }
        elif upperLowerEq == 'equality':
            maps = {'ALWAYS':self._ebndmap,
                    'AT_START':self._ebndmapStart,
                    'AT_END':self._ebndmapEnd,
                    }
        else:
            raise Exception('the "impossible" happened, upperLowerEq was not in '+\
                                "['upper','lower','equality']")

        def insertWithErr(bndmap):
           if name in bndmap:
               raise Exception(upperOrLower+' bound for '+name+' is already set to '+ \
                                   str(bndmap[name])+' but you tried to change it to '+str(bnd))
           bndmap[name] = bnd

        if when is None:
            insertWithErr(maps['ALWAYS'])
        elif when is "AT_END":
            insertWithErr(maps['AT_END'])
        elif when is "AT_START":
            insertWithErr(maps['AT_START'])
        else:
            raise Exception('the "impossible" happened: unrecognized "when": '+str(when))

    def constrain(self, *args, **kwargs):
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

        assert comparison in ['==','<=','>='], 'THE "IMPOSSIBLE" HAPPENED: comparison "'+str(comparison)+\
            '" is not "==", "<=", or ">="'

        # detect simple constraints in x and u
        def maybeAddBoxConstraint():
            inputs = C.veccat([self.dae.xVec(),self.dae.uVec()])
            # make sure only x and u are in rhs,lhs
            rml = rhs-lhs
            f = C.SXFunction([inputs],[rml])
            f.init()
            if len(f.getFree()) != 0:
                return

            # take jacobian of rhs-lhs
            jac = C.jacobian(rml,inputs)
            # fail if any jacobian element is not constant
            coeffs = {}
            for j in range(inputs.size()):
                if not jac[0,j].toScalar().isZero():
                    if not jac[0,j].toScalar().isConstant():
                        return
                    coeffs[j] = jac[0,j]
            if len(coeffs) == 0:
                raise Exception("constraint has no design variables in it")
            if len(coeffs) > 1:
                self.debug("found linear constraint that is not box constraint")
                return

            # alright, we've found a box constraint!
            j = coeffs.keys()[0]
            coeff = coeffs[j]
            name = (self.dae.xNames()+self.dae.uNames())[j]
            [f0] = f.eval([0*inputs])
            # if we just divided by a negative number (coeff), flip the comparison
            if not coeff.toScalar().isNonNegative():
                # lhs       `cmp`       rhs
                # 0         `cmp`       rhs - lhs
                # 0         `cmp`       coeff*x + f0
                # -f0       `cmp`       coeff*x
                # -f0/coeff `FLIP(cmp)` x
                if comparison == '>=':
                    newComparison = '<='
                elif comparison == '<=':
                    newComparison = '>='
                else:
                    newComparison = comparison
            else:
                newComparison = comparison

            c = -f0/coeff
            self.debug('found linear constraint: '+str(c)+' '+newComparison+' '+name)
            if newComparison == '==':
                self._bound( name, c, 'equality', when=when)
            elif newComparison == '<=':
                self._bound( name, c, 'lower', when=when)
            elif newComparison == '>=':
                self._bound( name, c, 'upper', when=when)
            else:
                raise Exception('the "impossible" happened, comparison "'+str(comparison)+
                                "\" not in ['==','>=','<=']")
            return 'found box constraint'

        if maybeAddBoxConstraint() == 'found box constraint':
            return

        if when is None:
            self._constraints.append((lhs,comparison,rhs))
        elif when is "AT_END":
            self._constraintsEnd.append((lhs,comparison,rhs))
        elif when is "AT_START":
            self._constraintsStart.append((lhs,comparison,rhs))
        else:
            raise Exception("\"when\" must be 'AT_START' or 'AT_END', leaving it blank means always, you put: "+str(when))

    def minimizeLsq(self, obj):
        if isinstance(obj, list):
            obj = C.veccat(obj)
        C.makeDense(obj)
        shape = obj.shape
        assert shape[0] == 1 or shape[1] == 1, 'objective cannot be matrix, got shape: '+str(shape)
        assert not hasattr(self, '_minLsq'), 'you can only call minimizeLsq once'
        self._minLsq = obj

    def minimizeLsqEndTerm(self, obj):
        if isinstance(obj, list):
            obj = C.veccat(obj)
        C.makeDense(obj)
        shape = obj.shape
        assert shape[0] == 1 or shape[1] == 1, 'objective cannot be matrix, got shape: '+str(shape)
        assert not hasattr(self, '_minLsqEndTerm'), 'you can only call minimizeLsqEndTerm once'
        self._minLsqEndTerm = obj

    def exportCode(self, ocpOptions, integratorOptions, codegenOptions, phase1Options):
        assert isinstance(ocpOptions, OcpExportOptions)
        assert isinstance(integratorOptions, RtIntegratorOptions)
        return exportOcp.exportOcp(self, ocpOptions, integratorOptions,
                                   codegenOptions, phase1Options)


class Mpc(Ocp):
    @property
    def yxNames(self):
        return self._yxNames
    @property
    def yuNames(self):
        return self._yuNames
    @property
    def yx(self):
        return self._yx
    @property
    def yu(self):
        return self._yu
    def __init__(self, dae, N=None, ts=None):
        Ocp.__init__(self, dae, N=N, ts=ts)
        self.hashPrefix = 'mpc'

        self._yxNames = dae.xNames()
        self._yuNames = dae.uNames()

        self._yx  = C.veccat( [self[n] for n in self.yxNames] )
        self._yu  = C.veccat( [self[n] for n in self.yuNames] )
        Ocp.minimizeLsq(self, C.veccat([self.yx,self.yu]))
        Ocp.minimizeLsqEndTerm(self, self.yx)

    def minimizeLsq(self, obj):
        raise Exception("hey, you don't know this is Ocp, the LSQ to be minimized is [X,U]")
    def minimizeLsqEndTerm(self, obj):
        raise Exception("hey, you don't know this is Ocp, the terminal cost is LSQ of X")


class Mhe(Ocp):
    @property
    def yxNames(self):
        return self._yxNames
    @property
    def yuNames(self):
        return self._yuNames
    @property
    def yx(self):
        return self._yx
    @property
    def yu(self):
        return self._yu
    def __init__(self, dae, N=None, ts=None, yxNames=None, yuNames=None):
        Ocp.__init__(self, dae, N=N, ts=ts)
        self.hashPrefix = 'mhe'

        if yxNames is None:
            yxNames = dae.xNames()
        if yuNames is None:
            yuNames = dae.uNames()
        if not isinstance(yxNames,list):
            raise Exception("If you decide to provide measurements, "+\
                            "you have to provide them as a list of strings")
        if not isinstance(yuNames,list):
            raise Exception("If you decide to provide end measurements, "+\
                            "you have to provide them as a list of strings")
        self._yxNames = yxNames
        self._yuNames = yuNames

        self._yx  = C.veccat( [self[n] for n in self.yxNames] )
        self._yu  = C.veccat( [self[n] for n in self.yuNames] )

        assert not C.dependsOn(self.yx, self.dae.uVec()), "error: x measurement depends on u"
        assert not C.dependsOn(self.yu, self.dae.xVec()), "error: u measurement depends on x"

        Ocp.minimizeLsq(self,C.veccat([self.yx, self.yu]))
        Ocp.minimizeLsqEndTerm(self,self.yx)

    def minimizeLsq(self, obj):
        raise Exception("hey, you don't know this is Ocp")
    def minimizeLsqEndTerm(self, obj):
        raise Exception("hey, you don't know this is Ocp")
