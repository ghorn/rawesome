import zmq
import time
import os
import numpy
import copy
from collections import Counter

import casadi as C

import kite_pb2
import model
import joy

#tc0 = 2*389.970797939731

x0 = C.DMatrix( [ 1.154244772411
                , -0.103540608242
                , -0.347959211327
                , 0.124930983341
                , 0.991534857363
                , 0.035367725910
                , 0.316039689643
                , -0.073559821379
                , 0.945889986864
                , 0.940484536806
                , -0.106993361072
                , -0.322554269411
                , 0.000000000000
                , 0.000000000000
                , 0.000000000000
                , 0.137035790811
                , 3.664945343102
                , -1.249768772258
                , 0.000000000000
                , 3.874600000000
                ])

r = 1.2

ts = 0.02

def toProto(x,u):
    cs = kite_pb2.CarouselState()

    cs.kiteXyz.x = x.at(0)
    cs.kiteXyz.y = x.at(1)
    cs.kiteXyz.z = x.at(2)

    cs.kiteDcm.r11 = x.at(3)
    cs.kiteDcm.r12 = x.at(4)
    cs.kiteDcm.r13 = x.at(5)

    cs.kiteDcm.r21 = x.at(6)
    cs.kiteDcm.r22 = x.at(7)
    cs.kiteDcm.r23 = x.at(8)

    cs.kiteDcm.r31 = x.at(9)
    cs.kiteDcm.r32 = x.at(10)
    cs.kiteDcm.r33 = x.at(11)

    cs.delta = x.at(18)
    cs.ddelta = x.at(19)
    
    cs.tc = u.at(0)
    cs.u0 = u.at(1)
    cs.u1 = u.at(2)
    return cs

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

        self.xBnd = [[None for n in range(0,self.nSteps)] for m in xNames]
        self.uBnd = [[None for n in range(0,self.nSteps)] for m in uNames]
        self.pBnd = [[None for n in range(0,self.nSteps)] for m in pNames]

        self.boundFcns = {}

        def makeBoundFcn(names,boundMat):
            for k,name in enumerate(names):
                self.boundFcns[name] = (k,boundMat)
        
        makeBoundFcn(xNames,self.xBnd)
        makeBoundFcn(uNames,self.uBnd)
        makeBoundFcn(pNames,self.pBnd)

    def bound(self,name,lbub,timestep=None):
        if timestep==None:
            for timestep in range(0,self.nSteps):
                self.bound(name,lbub,timestep)
            return
        
        if name not in self.boundFcns:
            raise ValueError("unrecognized variable name \""+name+"\"")
        #print "called bound("+name+","+str(lbub)+","+str(timestep)+")"
        (k,boundMat) = self.boundFcns[name]
        bnd0 = boundMat[k][timestep]
        if bnd0 != None:
            print "WARNING: bound for \""+name+"\" at timestep "+str(timestep)+" being changed from "+str(bnd0)+" to "+str(lbub)
        boundMat[k][timestep] = lbub

    def _concat(self,blah):
        import itertools
        chain = itertools.chain(*blah)
        return list(chain)

    def _concatBnds(self):
        return self._concat(self.xBnd)+self._concat(self.uBnd)+self._concat(self.pBnd)
    
    def getLb(self):
        return [x[0] for x in self._concatBnds()]

    def getUb(self):
        return [x[1] for x in self._concatBnds()]

def main():
    nSteps = 10
    endTime = C.ssym('endTime')

    print "creating model"
    (dae, others) = model.model((endTime,nSteps))
    dae.init()

#    # compile model code
#    print "generating model code"
#    t0 = time.time()
#    dae.generateCode("dae.c")
#    print "took "+str(time.time()-t0)+" seconds to generate code"
#    t0 = time.time()
#    os.system("gcc -fPIC -O2 -shared dae.c -o dae.so")
#    print "took "+str(time.time()-t0)+" seconds to compile code"
#    dae_ext = C.ExternalFunction("./dae.so")
#    dae_ext.init()
#    dae = dae_ext

    nStates = others['xVec'].size()
    nActions = others['uVec'].size()
    nParams = others['pVec'].size()

    assert(nStates==dae.inputSX(C.DAE_X).size())
    assert(nActions+nParams==dae.inputSX(C.DAE_P).size())
    
    # make the integrator
    print "creating integrator"
    integrator = C.IdasIntegrator(dae)
    integrator.setOption("reltol",1e-6)
    integrator.setOption("abstol",1e-8)
    integrator.setOption("t0",0)
    integrator.setOption("tf",1)
    integrator.setOption('name','integrator')
    integrator.init()

    # make the OCP
    print "setting up OCP"
    states = C.msym("x" ,nStates,nSteps)
    actions = C.msym("u",nActions,nSteps)
    params = C.msym("p",nParams)
    constraints = Constraints()

    constraints.addDynamicsConstraints(integrator, states, actions, params)

    # constrain invariants
    def invariantErrs():
        dcm = C.horzcat( [ C.veccat([others['xDict']['e11'], others['xDict']['e21'], others['xDict']['e31']])
                         , C.veccat([others['xDict']['e12'], others['xDict']['e22'], others['xDict']['e32']])
                         , C.veccat([others['xDict']['e13'], others['xDict']['e23'], others['xDict']['e33']])
                         ] ).trans()
        err = C.mul(dcm.trans(), dcm)
        dcmErr = C.veccat([ err[0,0]-1, err[1,1]-1, err[2,2]-1, err[0,1], err[0,2], err[1,2] ])
        f = C.SXFunction( [others['xVec'],others['uVec'],others['pVec']]
                        , [others['c'],others['cdot'],dcmErr]
                        )
        f.setOption('name','invariant errors')
        f.init()
        return f
    
    f = invariantErrs()
    [c0,cdot0,dcmError0] = f.call([states[:,0],actions[:,0],params])
    constraints.add(c0,'==',r)
    constraints.add(cdot0,'==',0)
    constraints.add(dcmError0,'==',0)

    # make it periodic
    constraints.add(states[:,0],'==',states[:,-1])
    constraints.add(actions[:,0],'==',actions[:,-1])

    # bounds
    bounds = Bounds(nSteps, others['xNames'], others['uNames'], others['pNames'])
    bounds.bound('u1',(-0.04,0.04))
    bounds.bound('u2',(-0.1,0.1))
    
    bounds.bound('y',(-3,3),3)
    bounds.bound('y',(-4,4))

    # make the solver
    designVars = C.veccat( [C.flatten(states), C.flatten(actions), C.flatten(params)] )
    
    # objective function
    obj = C.sumAll(actions*actions)
    f = C.MXFunction([designVars], [obj])

    # constraint function
    g = C.MXFunction([designVars], [constraints.getG()])

    # solver
    solver = C.IpoptSolver(f, g)
    solver.init()

    solver.setInput(constraints.getLb(), C.NLP_LBG)
    solver.setInput(constraints.getUb(), C.NLP_UBG)

    solver.setInput(bounds.getLb(), C.NLP_LBX)
    solver.setInput(bounds.getUb(), C.NLP_UBX)

    solver.setInput(numpy.zeros(designVars.size()), C.NLP_X_INIT)

    solver.solve()

if __name__=='__main__':
    main()
