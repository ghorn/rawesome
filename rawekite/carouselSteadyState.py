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
from rawe.ocputils import Constraints

pi = C.pi

def makeOrthonormal(g, R):
    g.add(C.mul(R[0, :], R[0, :].T), '==', 1, tag = ('R1[0]: e1^T * e1 == 1', None))
    g.add(C.mul(R[1, :], R[0, :].T), '==', 0, tag = ('R1[0]: e2^T * e1 == 0', None))
    g.add(C.mul(R[1, :], R[1, :].T), '==', 1, tag = ('R1[0]: e2^T * e2 == 1', None))
    rhon = C.cross(R[0, :], R[1, :]) - R[2, :]
    g.add(rhon[0], '==', 0, tag = ('R1[0]: ( e1^T X e2 - e3 )[0] == 0', None))
    g.add(rhon[2], '==', 0, tag = ('R1[0]: ( e1^T X e2 - e3 )[1] == 0', None))
    g.add(rhon[1], '==', 0, tag = ('R1[0]: ( e1^T X e2 - e3 )[2] == 0', None))

def getSteadyState(dae, conf, omega0, r0, ref_dict = {}):
    # make steady state model
    g = Constraints()
    g.add(dae.getResidual(), '==', 0, tag = ('dae residual', None))
    def constrainInvariantErrs():
        R_c2b = dae['R_c2b']
        makeOrthonormal(g, R_c2b)
        g.add(dae['c'], '==', 0, tag = ('c(0) == 0', None))
        g.add(dae['cdot'], '==', 0, tag = ('cdot( 0 ) == 0', None))
    constrainInvariantErrs()

    # Rotational velocity time derivative
    g.add(C.mul(dae['R_c2b'].T, dae['w_bn_b']) - C.veccat([0, 0, omega0]) , '==', 0, tag =
                       ("Rotational velocities", None))

    for name in ['alpha_deg', 'beta_deg', 'cL']:
        if name in ref_dict:
            g.addBnds(dae[name], ref_dict[name], tag = (name, None))

    dvs = C.veccat([dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), dae.xDotVec()])
    obj = 0
    for name in ['aileron', 'elevator', 'rudder', 'flaps']:
        if name in dae:
            obj += dae[name] ** 2
    ffcn = C.SXFunction([dvs], [obj])
    gfcn = C.SXFunction([dvs], [g.getG()])
    ffcn.init()
    gfcn.init()

    guess = {'x':r0, 'y':0, 'z':0,
             'r':r0, 'dr':0,
             'e11':0, 'e12':-1, 'e13':0,
             'e21':0, 'e22':0, 'e23':1,
             'e31':-1, 'e32':0, 'e33':0,
             'dx':0, 'dy':0, 'dz':0,
             'w_bn_b_x':0, 'w_bn_b_y':omega0, 'w_bn_b_z':0,
             'ddelta':omega0,
             'cos_delta':1, 'sin_delta':0,
             'aileron':0, 'elevator':0, 'rudder':0, 'flaps':0,
             'daileron':0, 'delevator':0, 'drudder':0, 'dflaps':0,
             'nu':100, 'motor_torque':0,
             'dmotor_torque':0, 'ddr':0,
             'dddr':0.0, 'w0':0.0}
    dotGuess = {'x':0, 'y':0, 'z':0, 'dx':0, 'dy':0, 'dz':0,
                'r':0, 'dr':0,
                'e11':0, 'e12':0, 'e13':0,
                'e21':0, 'e22':0, 'e23':0,
                'e31':0, 'e32':0, 'e33':0,
                'w_bn_b_x':0, 'w_bn_b_y':0, 'w_bn_b_z':0,
                'ddelta':0,
                'cos_delta':0, 'sin_delta':omega0,
                'aileron':0, 'elevator':0, 'rudder':0, 'flaps':0,
                'motor_torque':0, 'ddr':0}

    guessVec = C.DMatrix([guess[n] for n in dae.xNames() + dae.zNames() + dae.uNames() + dae.pNames()] +
                         [dotGuess[n] for n in dae.xNames()])

    bounds = {'x':(-0.1 * r0, r0 * 2), 'y':(-1.1 * r0, 1.1 * r0), 'z':(-1.1 * r0, 1.1 * r0),
             'dx':(0, 0), 'dy':(0, 0), 'dz':(0, 0),
             'r':(r0, r0), 'dr':(-100, 100),
             'e11':(-2, 2), 'e12':(-2, 2), 'e13':(-2, 2),
             'e21':(-2, 2), 'e22':(-2, 2), 'e23':(-2, 2),
             'e31':(-2, 2), 'e32':(-2, 2), 'e33':(-2, 2),
             'w_bn_b_x':(-50, 50), 'w_bn_b_y':(-50, 50), 'w_bn_b_z':(-50, 50),
             'ddelta':(omega0, omega0),
             'cos_delta':(1, 1), 'sin_delta':(0, 0),
             'aileron':(-0.2, 0.2), 'elevator':(-0.2, 0.2), 'rudder':(-0.2, 0.2), 'flaps':(-0.2, 0.2),
             'daileron':(0, 0), 'delevator':(0, 0), 'drudder':(0, 0), 'dflaps':(0, 0),
             'nu':(0, 3000), 'motor_torque':(-1000, 1000),
             'ddr':(0, 0),
             'dmotor_torque':(0, 0), 'dddr':(0, 0), 'w0':(0, 0)}

    if ref_dict is not None:
        for name in ref_dict: bounds[name] = ref_dict[name]

    print bounds

    dotBounds = {'x':(-1, 1), 'y':(-1, 1), 'z':(-1, 1),
                 'dx':(0, 0), 'dy':(0, 0), 'dz':(0, 0),
                 'r':(-100, 100), 'dr':(-1, 1),
                 'e11':(-50, 50), 'e12':(-50, 50), 'e13':(-50, 50),
                 'e21':(-50, 50), 'e22':(-50, 50), 'e23':(-50, 50),
                 'e31':(-50, 50), 'e32':(-50, 50), 'e33':(-50, 50),
                 'w_bn_b_x':(0, 0), 'w_bn_b_y':(0, 0), 'w_bn_b_z':(0, 0),
                 'ddelta':(0, 0),
                 'cos_delta':(-1, 1), 'sin_delta':(omega0 - 1, omega0 + 1),
                 'aileron':(-1, 1), 'elevator':(-1, 1), 'rudder':(-1, 1), 'flaps':(-1, 1),
                 'motor_torque':(-1000, 1000), 'ddr':(-100, 100)}
    boundsVec = [bounds[n] for n in dae.xNames() + dae.zNames() + dae.uNames() + dae.pNames()] + \
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


    solver = C.IpoptSolver(ffcn, gfcn)
#    def addCallback():
#        nd = len(boundsVec)
#        nc = g.getLb().size()
#        c = C.PyFunction( MyCallback(), C.nlpsolverOut(x=C.sp_dense(nd,1), f=C.sp_dense(1,1), lam_x=C.sp_dense(nd,1), lam_p = C.sp_dense(0,1), lam_g = C.sp_dense(nc,1), g = C.sp_dense(nc,1) ), [C.sp_dense(1,1)] )
#        c.init()
#        solver.setOption("iteration_callback", c)
#    addCallback()
    solver.setOption('max_iter', 10000)
#    solver.setOption('tol',1e-14)
#    solver.setOption('suppress_all_output','yes')
#    solver.setOption('print_time',False)
    solver.init()

    solver.setInput(g.getLb(), 'lbg')
    solver.setInput(g.getUb(), 'ubg')
    solver.setInput(guessVec, 'x0')
    lb, ub = zip(*boundsVec)
    solver.setInput(C.DMatrix(lb), 'lbx')
    solver.setInput(C.DMatrix(ub), 'ubx')

    solver.solve()
    ret = solver.getStat('return_status')
    assert ret in ['Solve_Succeeded', 'Solved_To_Acceptable_Level'], 'Solver failed: ' + ret

#    publisher.close()
#    context.destroy()
    xOpt = solver.output('x')
    k = 0
    sol = {}
    for name in dae.xNames() + dae.zNames() + dae.uNames() + dae.pNames():
        sol[name] = xOpt[k].at(0)
        k += 1
#        print name+':\t',sol[name]
    dotSol = {}
    for name in dae.xNames():
        dotSol[name] = xOpt[k].at(0)
        k += 1
#        print 'DDT('+name+'):\t',dotSol[name]
    return sol, dotSol
