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

import matplotlib.pyplot as plt
import numpy

import rawe
import casadi as C

# specify the dae
def makeDae():
    dae = rawe.dae.Dae()

    [pos,vel] = dae.addX( ["pos","vel"] )
    force = dae.addU( "force" )

    dae.setResidual([dae.ddt('pos') - vel,
                     dae.ddt('vel') - (force - 0.14*pos - 0*0.2*vel)])
    return dae

def makeMpc(dae, N, ts):
    mpc = rawe.Ocp(dae, N=N, ts=ts)
    mpc.constrain(-2.5, '<=', mpc['force'], '<=', 2.5)

    mpc.minimizeLsq([mpc['pos'],mpc['vel'],mpc['force']])
    mpc.minimizeLsqEndTerm([mpc['pos'],mpc['vel']])

#    cgOpts = {'CXX':'clang++', 'CC':'clang'}
    cgOpts = {'CXX':'g++', 'CC':'gcc'}
#    cgOpts = {'CXX':'icpc', 'CC':'icc'}
    intOpts = rawe.RtIntegratorOptions()
    intOpts['INTEGRATOR_TYPE'] = 'INT_IRK_GL4'
    intOpts['NUM_INTEGRATOR_STEPS'] = 5
    intOpts['LINEAR_ALGEBRA_SOLVER'] = 'GAUSS_LU'

    ocpOpts = rawe.OcpExportOptions()
    ocpOpts['HESSIAN_APPROXIMATION'] = 'GAUSS_NEWTON'
    ocpOpts['DISCRETIZATION_TYPE'] = 'MULTIPLE_SHOOTING'
    ocpOpts['QP_SOLVER'] = 'QP_QPOASES'
    ocpOpts['SPARSE_QP_SOLUTION'] = 'CONDENSING'
#    ocpOpts['SPARSE_QP_SOLUTION'] = 'FULL_CONDENSING'
#    ocpOpts['QP_SOLVER'] = 'QP_FORCES'
#    ocpOpts['SPARSE_QP_SOLUTION'] = 'SPARSE_SOLVER'
    ocpOpts['FIX_INITIAL_STATE'] =      True
    ocpOpts['HOTSTART_QP'] =            True
    ocpOpts['GENERATE_MATLAB_INTERFACE'] = True
    return rawe.OcpRT(mpc, ocpOptions=ocpOpts, integratorOptions=intOpts,
                       codegenOptions=cgOpts)

def makeMhe(dae, N, ts):
    mhe = rawe.ocp.Ocp(dae, N=N, ts=ts)

    mhe.minimizeLsq([mhe['pos'],mhe['vel']])
    mhe.minimizeLsqEndTerm([mhe['pos'],mhe['vel']])

#    cgOpts = {'CXX':'clang++', 'CC':'clang'}
    cgOpts = {'CXX':'g++', 'CC':'gcc'}
#    cgOpts = {'CXX':'icpc', 'CC':'icc'}
    return rawe.OcpRT(mhe, codegenOptions=cgOpts)


if __name__=='__main__':
    N = 70
    ts = 0.2

    A = 0.5
    omega = 0.3

    dae = makeDae()
    mpc = makeMpc(dae, N, ts)

    sim = rawe.RtIntegrator(dae, ts=ts)
#    sim = rawe.sim.Sim(dae, ts=ts)

    print '='*80

    # set the mpc weights
    xRms = 0.1
    vRms = 1.5
    fRms = 5.0
    mpc.S[0,0] = (1.0/xRms)**2/N
    mpc.S[1,1] = (1.0/vRms)**2/N
    mpc.S[2,2] = (1.0/fRms)**2/N
    mpc.SN[0,0] = (1.0/0.01)**2
    mpc.SN[1,1] = (1.0/0.01)**2

    # initial guess
    tk = 0
    for k in range(N+1):
        mpc.x[k,0] = A*C.sin(omega*tk)
        mpc.x[k,1] = A*C.cos(omega*tk)*omega
        tk += ts

    # set tracking trajectory
    t = 0
    refLog = []
    def setTrajectory(ocp, t0):
        tk = t0
        for k in range(N):
            ocp.y[k,0] = A*C.sin(omega*tk)
            ocp.y[k,1] = A*C.cos(omega*tk)*omega
            tk += ts
        ocp.yN[0] = A*C.sin(omega*tk)
        ocp.yN[1] = A*C.cos(omega*tk)*omega

        refLog.append(numpy.copy(numpy.vstack((ocp.y[:,:2], ocp.yN.T))))
    setTrajectory(mpc, t)

    # run a sim
    x = {'pos':mpc.x[0,0],'vel':mpc.x[0,1]}
    u = {'force':0}

    ytraj = []
    xtraj = []
    xestTraj = []
    utraj = []
    kkts = []
    objs = []
    t = 0
    for k in range(200):
        mpc.x0[0] = x['pos'] # will be mhe output
        mpc.x0[1] = x['vel'] # will be mhe output

        # command mpc to follow sin wave
        setTrajectory(mpc, t)

        # log stuff
        xtraj.append((x['pos'], x['vel']))
        ytraj.append((mpc.y[0,0], mpc.y[0,1]))
        xestTraj.append((mpc.x0[0], mpc.x0[1]))

        # run mpc
        for j in range(2):
            mpc.preparationStep()
            mpc.feedbackStep()
            print "timestep",k,"sqp iteration",j,"\tkkts:",mpc.getKKT(),\
                "\tobjective:",mpc.getObjective(),\
                "\tpreparation time:",mpc.preparationTime,"\tfeedback time:",mpc.feedbackTime


        u['force'] = mpc.u[0,0]

        # log
        utraj.append(u['force'])
        kkts.append(mpc.getKKT() + 1e-100)
        objs.append(mpc.getObjective() + 1e-100)

        # simulate
        x = sim.step(x, u, {})
        t += ts

    # plot results
    plt.figure()
    plt.subplot(311)
    plt.plot([xt[0] for xt in ytraj])
    plt.plot([xt[0] for xt in xtraj])
    plt.legend(['x command','x actual'])

    plt.subplot(312)
    plt.plot([xt[1] for xt in ytraj])
    plt.plot([xt[1] for xt in xtraj])
    plt.legend(['v command','v actual'])

    plt.subplot(313)
    plt.plot(utraj)
    plt.legend(['u'])


    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for k,y in enumerate(refLog):
        n = y[:,0].size
        xs = range(k,k+n)
        ax.plot(xs, k*numpy.ones(n), zs=y[:,0] )
    ax.set_xlabel("time")
    ax.set_ylabel("iteration")
    ax.set_zlabel("reference")
    plt.title('pos ref')

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for k,y in enumerate(refLog):
        n = y[:,1].size
        xs = range(k,k+n)
        ax.plot(xs, k*numpy.ones(n), zs=y[:,1] )
    ax.set_xlabel("time")
    ax.set_ylabel("iteration")
    ax.set_zlabel("reference")
    plt.title('vel ref')

    plt.figure()
    plt.subplot(211)
    plt.semilogy(kkts)
    plt.title('kkts')
    plt.subplot(212)
    plt.semilogy(objs)
    plt.title('objectives')

    plt.show()
