#! /usr/bin/env ipython

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

if __name__=='__main__':
    # specify the dae
    dae = rawe.dae.Dae()

    [pos,vel] = dae.addX( ["pos","vel"] )
    force = dae.addU( "force" )
    springK = dae.addP( "springK" )

    dae.setResidual([dae.ddt('pos') - vel,
                     dae.ddt('vel') - (force - springK*pos - 0.2*vel)])

    # specify the Ocp
    N = 100
    mpc = rawe.Ocp(dae, N=N, ts=0.2, yxNames=['pos', 'vel'], yuNames=['force'])

#    mpc.constrain(mpc['pos'], '==', 0, when='AT_START')
#    mpc.constrain(mpc['vel'], '==', 0, when='AT_START')

    mpc.constrain(mpc['pos'], '==', 0.1, when='AT_END')
    mpc.constrain(mpc['vel'], '==', 0, when='AT_END')

    mpc.constrain(mpc['vel'], '<=', 0.2)
    mpc.constrain(-0.1, '<=', mpc['force'], '<=', 0.1)

#    cgOpts = {'CXX':'clang++', 'CC':'clang'}
    cgOpts = {'CXX':'g++', 'CC':'gcc'}
#    cgOpts = {'CXX':'icpc', 'CC':'icc'}
    phase1Opts = {'CXX':'g++'}
    intOpts = rawe.RtIntegratorOptions()
    intOpts['INTEGRATOR_TYPE'] = 'INT_IRK_RIIA3'
    intOpts['NUM_INTEGRATOR_STEPS'] = 5

    ocpOpts = rawe.OcpExportOptions()
    ocpOpts['QP_SOLVER'] = "QP_QPOASES"
    ocpOpts['SPARSE_QP_SOLUTION'] = "CONDENSING"
#    ocpOpts['SPARSE_QP_SOLUTION'] = "FULL_CONDENSING"
#    ocpOpts['SPARSE_QP_SOLUTION'] = "SPARSE_SOLVER"
#    ocpOpts['QP_SOLVER'] = "QP_QPDUNES"
    ocpOpts["FIX_INITIAL_STATE"] = True
    ocpOpts["HOTSTART_QP"] = True
#   ocpOpts['GENERATE_MATLAB_INTERFACE'] = True

    ocpRt = rawe.OcpRT(mpc, ocpOptions=ocpOpts, integratorOptions=intOpts,
                       codegenOptions=cgOpts, phase1Options=phase1Opts)
    print '='*80

    # set the cost hessians
#    xRms = 0.2
#    vRms = 0.2
#    fRms = 20
    ocpRt.S[0,0] = 1.0#/(xRms*N)**2
    ocpRt.S[1,1] = 1.0#/(vRms*N)**2
    ocpRt.S[2,2] = 1.0#/(fRms*N)**2

    ocpRt.SN[0,0] = 1.0#/(xRms*N)**2
    ocpRt.SN[1,1] = 1.0#/(vRms*N)**2

    # make an initial guess
    ocpRt.x[0,0] = 0.0
    ocpRt.x[0,1] = 0.0
    for k in range(N):
        ocpRt.u[k] = 0.0
    ocpRt.p[:,0] = 3.0
    ocpRt.initializeNodesByForwardSimulation()

#    import pickle
#    (ocpRt.x, ocpRt.u) = pickle.load( open('rt_initial_guess.dat','rb'))

    # plot initial guess
    plt.figure()
    plt.subplot(331)
    plt.plot(ocpRt.x[:,0])
    plt.title('x initial guess')
    plt.subplot(334)
    plt.plot(ocpRt.x[:,1])
    plt.title('v initial guess')
    plt.subplot(337)
    plt.plot(ocpRt.u[:,0])
    plt.title('u initial guess')

    # set tracking trajectory
    ocpRt.yN[0] = 0
    ocpRt.yN[1] = 0

    # iterate
    kkts = []
    objs = []

    xHistory = [numpy.array(ocpRt.x)]
    uHistory = [numpy.array(ocpRt.u)]
    for k in range(4):
        ocpRt.preparationStep()
        ocpRt.feedbackStep()
        print "sqp iteration",k,"\tkkts:",ocpRt.getKKT(),"\tobjective:",ocpRt.getObjective(),\
            "\tpreparation time:",ocpRt.preparationTime,"\tfeedback time:",ocpRt.feedbackTime
        xHistory.append(numpy.array(ocpRt.x))
        uHistory.append(numpy.array(ocpRt.u))

        kkts.append(ocpRt.getKKT() + 1e-200)
        objs.append(ocpRt.getObjective() + 1e-200)

    ocpRt.writeStateTxtFiles()

#    print "1:"
#    print xHistory[1].T - xHistory[2].T
#    print uHistory[1].T - uHistory[2].T
#    import sys; sys.exit()

#    import pickle
#    print "saving to rt_initial_guess.dat"
#    pickle.dump( (ocpRt.x, ocpRt.u), open('rt_initial_guess.dat','wb'))

    # plot results
    plt.subplot(332)
    plt.plot(ocpRt.x[:,0])
    plt.title('x after feedback')
    plt.subplot(335)
    plt.plot(ocpRt.x[:,1])
    plt.title('v after feedback')
    plt.subplot(338)
    plt.plot(ocpRt.u[:,0])
    plt.title('u after feedback')

    plt.subplot(233)
    plt.semilogy(kkts)
    plt.title('kkt convergence')
    plt.subplot(236)
    plt.semilogy(objs)
    plt.title('objective convergence')
    plt.show()
