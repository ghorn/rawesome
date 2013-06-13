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

import zmq
import casadi as C
import numpy

from rawe.ocputils import Constraints

pi = C.pi

def makeOrthonormal(g,R):
         g.add(C.mul(R[0,:],R[0,:].T),'==',1,  tag=('R1[0]: e1^T * e1 == 1',None))
         g.add(C.mul(R[1,:],R[0,:].T),'==',0,  tag=('R1[0]: e2^T * e1 == 0',None))
         g.add(C.mul(R[1,:],R[1,:].T),'==',1,  tag=('R1[0]: e2^T * e2 == 1',None))
         rhon = C.cross(R[0,:],R[1,:]) - R[2,:]
         g.add(rhon[0],'==',0,  tag=('R1[0]: ( e1^T X e2 - e3 )[0] == 0',None))
         g.add(rhon[2],'==',0,  tag=('R1[0]: ( e1^T X e2 - e3 )[1] == 0',None))
         g.add(rhon[1],'==',0,  tag=('R1[0]: ( e1^T X e2 - e3 )[2] == 0',None))

def getSteadyState(dae,r0,v0):
    # make steady state model
    g = Constraints()
    g.add(dae.getResidual(),'==',0,tag=('dae residual',None))
    def constrainInvariantErrs():
        R_n2b = dae['R_n2b']
        makeOrthonormal(g, R_n2b)
        g.add(dae['c'], '==', 0, tag=('c(0)==0',None))
        g.add(dae['cdot'], '==', 0, tag=('cdot(0)==0',None))
    constrainInvariantErrs()

    # constrain airspeed
    g.add(dae['airspeed'], '>=', v0, tag=('airspeed fixed',None))
    g.addBnds(dae['alpha_deg'], (4,10), tag=('alpha',None))
    g.addBnds(dae['beta_deg'], (-10,10), tag=('beta',None))
    
    dvs = C.veccat([dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), dae.xDotVec()])
#    ffcn = C.SXFunction([dvs],[sum([dae[n]**2 for n in ['aileron','elevator','y','z']])])
    obj = 0
    obj += (dae['cL']-0.5)**2
    obj += dae.ddt('w_bn_b_x')**2
    obj += dae.ddt('w_bn_b_y')**2
    obj += dae.ddt('w_bn_b_z')**2
    ffcn = C.SXFunction([dvs],[obj])
    gfcn = C.SXFunction([dvs],[g.getG()])
    ffcn.init()
    gfcn.init()

    guess = {'x':r0,'y':0,'z':-1,
             'r':r0,'dr':0,
             'e11':0, 'e12':-1, 'e13':0,
             'e21':0, 'e22':0, 'e23':1,
             'e31':-1, 'e32':0, 'e33':0,
             'dx':0,'dy':-20,'dz':0,
             'w_bn_b_x':0,'w_bn_b_y':0,'w_bn_b_z':0,
             'aileron':0,'elevator':0,'rudder':0,
             'daileron':0,'delevator':0,'drudder':0,
             'nu':300,'motor_torque':10,
             'dmotor_torque':0,'ddr':0,
             'dddr':0.0,'w0':10.0}
    dotGuess = {'x':0,'y':-20,'z':0,'dx':0,'dy':0,'dz':0,
                'r':0,'dr':0,
                'e11':0,'e12':0,'e13':0,
                'e21':0,'e22':0,'e23':0,
                'e31':0,'e32':0,'e33':0,
                'w_bn_b_x':0,'w_bn_b_y':0,'w_bn_b_z':0,
                'aileron':0,'elevator':0,'rudder':0,
                'ddr':0}

    guessVec = C.DMatrix([guess[n] for n in dae.xNames()+dae.zNames()+dae.uNames()+dae.pNames()]+
                         [dotGuess[n] for n in dae.xNames()])

    bounds = {'x':(0.01,r0*2),'y':(0,0),'z':(-r0*0.2,-r0*0.2),
              'dx':(0,0),'dy':(-50,0),'dz':(0,0),
              'r':(r0,r0),'dr':(0,0),'ddr':(0,0),
              'e11':(-0.5,0.5),'e12':(-1.5,-0.5),'e13':(-0.5,0.5),
              'e21':(-0.5,0.5),'e22':(-0.5,0.5),'e23':(0.5,1.5),
              'e31':(-1.5,-0.5),'e32':(-0.5,0.5),'e33':(-0.5,0.5),
              'w_bn_b_x':(0,0),'w_bn_b_y':(0,0),'w_bn_b_z':(0,0),
#              'aileron':(-0.2,0.2),'elevator':(-0.2,0.2),'rudder':(-0.2,0.2),
              'aileron':(0,0),'elevator':(0,0),'rudder':(0,0),
              'daileron':(0,0),'delevator':(0,0),'drudder':(0,0),
              'nu':(0,3000),
              'dddr':(0,0),'w0':(10,10)}
    dotBounds = {'dz':(-C.inf,0)}#'dx':(-500,-500),'dy':(-500,500),'dz':(-500,500),
#                 'w_bn_b_x':(0,0),'w_bn_b_y':(0,0),'w_bn_b_z':(0,0),

    for name in dae.xNames():
        if name not in dotBounds:
            dotBounds[name] = (-C.inf, C.inf)
    boundsVec = [bounds[n] for n in dae.xNames()+dae.zNames()+dae.uNames()+dae.pNames()]+ \
                [dotBounds[n] for n in dae.xNames()]
    
#    gfcn.setInput(guessVec)
#    gfcn.evaluate()
#    ret = gfcn.output()
#    for k,v in enumerate(ret):
#        if math.isnan(v):
#            print 'index ',k,' is nan: ',g._tags[k]
#    import sys; sys.exit()

#    context   = zmq.Context(1)
#    publisher = context.socket(zmq.PUB)
#    publisher.bind("tcp://*:5563")
#    class MyCallback:
#        def __init__(self):
#            import rawekite.kiteproto as kiteproto
#            import rawekite.kite_pb2 as kite_pb2
#            self.kiteproto = kiteproto
#            self.kite_pb2 = kite_pb2
#            self.iter = 0
#        def __call__(self,f,*args):
#            x = C.DMatrix(f.input('x'))
#            sol = {}
#            for k,name in enumerate(dae.xNames()+dae.zNames()+dae.uNames()+dae.pNames()):
#                sol[name] = x[k].at(0)
#            lookup = lambda name: sol[name]
#            kp = self.kiteproto.toKiteProto(lookup,
#                                            lineAlpha=0.2)
#            mc = self.kite_pb2.MultiCarousel()
#            mc.horizon.extend([kp])
#            mc.messages.append("iter: "+str(self.iter))
#            self.iter += 1
#            publisher.send_multipart(["multi-carousel", mc.SerializeToString()])

    
    solver = C.IpoptSolver(ffcn,gfcn)
    def addCallback():
        nd = len(boundsVec)
        nc = g.getLb().size()
        c = C.PyFunction( MyCallback(), C.nlpsolverOut(x=C.sp_dense(nd,1), f=C.sp_dense(1,1), lam_x=C.sp_dense(nd,1), lam_p = C.sp_dense(0,1), lam_g = C.sp_dense(nc,1), g = C.sp_dense(nc,1) ), [C.sp_dense(1,1)] )
        c.init()
        solver.setOption("iteration_callback", c)
#    addCallback()
    solver.setOption('max_iter',10000)
    solver.setOption('expand',True)
#    solver.setOption('tol',1e-14)
#    solver.setOption('suppress_all_output','yes')
#    solver.setOption('print_time',False)
    solver.init()

    solver.setInput(g.getLb(),'lbg')
    solver.setInput(g.getUb(),'ubg')
    #guessVec = numpy.load('steady_state_guess.npy')
    solver.setInput(guessVec,'x0')
    lb,ub = zip(*boundsVec)
    solver.setInput(C.DMatrix(lb), 'lbx')
    solver.setInput(C.DMatrix(ub), 'ubx')

    solver.solve()
    ret = solver.getStat('return_status')
    assert ret in ['Solve_Succeeded','Solved_To_Acceptable_Level'], 'Solver failed: '+ret

#    publisher.close()
#    context.destroy()
    xOpt = solver.output('x')
    #numpy.save('steady_state_guess',numpy.squeeze(numpy.array(xOpt)))
    k = 0
    sol = {}
    for name in dae.xNames()+dae.zNames()+dae.uNames()+dae.pNames():
        sol[name] = xOpt[k].at(0)
        k += 1
#        print name+':\t',sol[name]
    dotSol = {}
    for name in dae.xNames():
        dotSol[name] = xOpt[k].at(0)
        k += 1
#        print 'DDT('+name+'):\t',dotSol[name]
    return sol, dotSol

