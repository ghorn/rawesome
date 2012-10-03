import zmq
#import time
#import os

import numpy
from numpy import pi
import copy

import casadi as C

import kite_pb2
import kiteproto
import ocp 
import model

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
x0=C.veccat([x0,C.sqrt(C.sumAll(x0[0:2]*x0[0:2])),0])

def main():
    nSteps = 15
    endTime = C.ssym('endTime')

    print "creating model"
    (dae, others, outputs) = model.model(-0.01,(endTime,nSteps))
    dae.init()

    nStates = others['xVec'].size()
    nActions = others['uVec'].size()
    nParams = others['pVec'].size()

    assert(nStates==dae.inputSX(C.DAE_X).size())
    assert(nActions+nParams==dae.inputSX(C.DAE_P).size())
    
    # make the integrator
    print "creating integrator"
    integrator = C.IdasIntegrator(dae)
    integrator.setOption("reltol",1e-7)
    integrator.setOption("abstol",1e-9)
    integrator.setOption("t0",0)
    integrator.setOption("tf",1)
    integrator.setOption('name','integrator')
    integrator.setOption("linear_solver_creator",C.CSparse)
#    integrator.setOption("linear_solver_creator",C.LapackLUDense)
    integrator.setOption("linear_solver","user_defined")
    integrator.init()

    # make the OCP
    print "setting up OCP"
    states = C.msym("x" ,nStates,nSteps)
    actions = C.msym("u",nActions,nSteps)
    params = C.msym("p",nParams)
    constraints = ocp.Constraints()

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
    
    [c0,cdot0,dcmError0] = invariantErrs().call([states[:,0],actions[:,0],params])
    constraints.add(c0,'==',0)
    constraints.add(cdot0,'==',0)
    constraints.add(dcmError0,'==',0)

    # make it periodic
    constraints.add(states[:18,0],'==',states[:18,-1])
    constraints.add(states[19,0],'==',states[19,-1])
    constraints.add(actions[:,0],'==',actions[:,-1])

    # bounds
    bounds = ocp.Bounds(others['xNames'], others['uNames'], others['pNames'], nSteps)
    bounds.setBound('aileron',(-0.04,0.04))
    bounds.setBound('elevator',(-0.1,0.1))
    
    bounds.setBound('x',(0,4))
    bounds.setBound('y',(-3,3))
    bounds.setBound('z',(-3,3))
    bounds.setBound('r',(1,2))
    bounds.setBound('dr',(-1,1))
    bounds.setBound('ddr',(0,0))
    bounds.setBound('r',(1.2,1.2),timestep=0)

    for e in ['e11','e21','e31','e12','e22','e32','e13','e23','e33']:
        bounds.setBound(e,(-1.1,1.1))

    for d in ['dx','dy','dz']:
        bounds.setBound(d,(-50,50))

    for w in ['w1','w2','w3']:
        bounds.setBound(w,(-4*pi,4*pi))

    bounds.setBound('delta',(-0.01,1.01*2*pi))
    bounds.setBound('ddelta',(-pi/4,8*pi))
    bounds.setBound('tc',(-200,600))
    bounds.setBound('tc',(389.970797939731,389.970797939731))
#    bounds.setBound('endTime',(0.5,2.0))
    bounds.setBound('endTime',(1.6336935276077966,1.6336935276077966))
    bounds.setBound('w0',(0,0))

    # boundary conditions
    bounds.setBound('delta',(0,0),timestep=0)
    bounds.setBound('delta',(2*pi,2*pi),timestep=nSteps-1)

    # make the solver
    designVars = C.veccat( [C.flatten(states), C.flatten(actions), C.flatten(params)] )
    
    # objective function
    dvs = ocp.DesignVars((others['xNames'],states), (others['uNames'],actions), (others['pNames'],params), nSteps)
    tc0 = 390
    obj = (C.sumAll(actions[0:2,:]*actions[0:2,:]) + 1e-10*C.sumAll((actions[2,:]-tc0)*(actions[2,:]-tc0)))*dvs.lookup('endTime')
    f = C.MXFunction([designVars], [obj])
    f.init()

    # constraint function
    g = C.MXFunction([designVars], [constraints.getG()])
    g.setOption('numeric_jacobian',False)
    g.init()

    context   = zmq.Context(1)
    publisher = context.socket(zmq.PUB)
    publisher.bind("tcp://*:5563")

    # callback function
    class MyCallback:
      def __init__(self):
        self.iter = 0 
      def __call__(self,f,*args):
          self.iter = self.iter + 1
          xOpt = f.input(C.NLP_X_OPT)

          xs,us,p = bounds.getTimestepsFromDvs(xOpt)
          kiteProtos = [kiteproto.toKiteProto(xs[k],us[k],p) for k in range(0,nSteps)]

          ko = kite_pb2.KiteOpt()
          ko.css.extend(list(kiteProtos))

          xup = bounds.devectorize(xOpt)
          ko.endTime = xup['endTime']
          ko.wind_x = xup['w0']
          ko.iters = self.iter
          publisher.send_multipart(["carousel-opt", ko.SerializeToString()])
    
    def makeCallback():
        nd = designVars.size()
        nc = constraints.getG().size()
        c = C.PyFunction( MyCallback(), C.nlpsolverOut(x_opt=C.sp_dense(nd,1), cost=C.sp_dense(1,1), lambda_x=C.sp_dense(nd,1), lambda_g = C.sp_dense(nc,1), g = C.sp_dense(nc,1) ), [C.sp_dense(1,1)] )
        c.init()
        return c
        

    # solver
    solver = C.IpoptSolver(f, g)
    solver.setOption("iteration_callback",makeCallback())
    solver.setOption("linear_solver","ma57")
    solver.init()

    solver.setInput(constraints.getLb(), C.NLP_LBG)
    solver.setInput(constraints.getUb(), C.NLP_UBG)

    lb,ub = bounds.get()
    solver.setInput(lb, C.NLP_LBX)
    solver.setInput(ub, C.NLP_UBX)

    # initial conditions
    guess = ocp.InitialGuess(others['xNames'], others['uNames'], others['pNames'], nSteps)
    guess.setXVec(x0)
    for k in range(0,nSteps):
        val = 2*pi*k/(nSteps-1)
        guess.setGuess('delta',val,timestep=k,quiet=True)

    guess.setGuess('aileron',0)
    guess.setGuess('elevator',0)
    guess.setGuess('tc',389.970797939731)
    guess.setGuess('endTime',1.6336935276077966)
    
    guess.setGuess('ddr',0)
    guess.setGuess('w0',0)


    solver.setInput(guess.vectorize(), C.NLP_X_INIT)

    solver.solve()

if __name__=='__main__':
    main()
