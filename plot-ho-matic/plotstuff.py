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
from loadprotos import loadProtos

protos = loadProtos("mhe-mpc-horizons_log.dat", start=1925, end=1975)

plt.figure()
plt.plot( [p.currentState.controls.dup for p in protos] )
plt.plot( [p.currentState.controls.dur for p in protos] )
plt.plot( [p.currentState.diffStates.x for p in protos] )
plt.plot( [p.currentState.diffStates.y for p in protos] )
plt.plot( [p.currentState.diffStates.z for p in protos] )
plt.legend(['dup','dur','x','y','z'])


for k in range(16,20):
    plt.figure()
    plt.plot( [h.dae.controls.dup for h in protos[k].referenceTrajectory] )
    plt.plot( [h.dae.controls.dur for h in protos[k].referenceTrajectory] )
    plt.plot( [h.dae.diffStates.x for h in protos[k].referenceTrajectory] )
    plt.plot( [h.dae.diffStates.y for h in protos[k].referenceTrajectory] )
    plt.plot( [h.dae.diffStates.z for h in protos[k].referenceTrajectory] )
    plt.legend(['dur','dup','x','y','z'])
    plt.title('reference traj #'+str(k))


#for k in range(16,20):
#    plt.figure()
#    plt.plot( [h.dae.controls.dup for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.controls.dur for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.x for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.y for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.z for h in protos[k].mpcHorizon] )
#    plt.legend(['dur','dup','x','y','z'])
#    plt.title('mpc horizon #'+str(k))


#for k in range(16,20):
#    plt.figure()
#    plt.plot( [h.dae.controls.dup for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.controls.dur for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.x for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.y for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.z for h in protos[k].mpcHorizon] )
#    plt.legend(['dur','dup','x','y','z'])
#    plt.title('mpc horizon #'+str(k))


#for k in range(16,20):
#    plt.figure()
#    plt.plot( [h.dae.diffStates.x      for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.y      for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.z      for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.dx     for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.dy     for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.dz     for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.e11    for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.e12    for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.e13    for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.e21    for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.e22    for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.e23    for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.e31    for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.e32    for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.e33    for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.wx     for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.wy     for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.wz     for h in protos[k].mpcHorizon] )
##    plt.plot( [h.dae.diffStates.delta  for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.ddelta for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.ur     for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.diffStates.up     for h in protos[k].mpcHorizon] )
#    plt.title('mpc horizon #'+str(k)+' diff states')


#for k in range(16,20):
#    plt.figure()
#    plt.plot( [h.dae.controls.dddelta for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.controls.dur for h in protos[k].mpcHorizon] )
#    plt.plot( [h.dae.controls.dup for h in protos[k].mpcHorizon] )
#    plt.legend(['dddelta','dur','dup'])
#    plt.title('mpc horizon #'+str(k)+' controls')
#
#
#for k in range(16,20):
#    plt.figure()
#    plt.plot( [h.measurementsX.uvC1M1_0 for h in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC1M1_1 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC1M2_0 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC1M2_1 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC1M3_0 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC1M3_1 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC2M1_0 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC2M1_1 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC2M2_0 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC2M2_1 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC2M3_0 for p in protos[k].measurementsHorizon] )
#    plt.plot( [h.measurementsX.uvC2M3_1 for p in protos[k].measurementsHorizon] )
#    plt.title('mpc horizon #'+str(k)+' camera meas')


################## measurements ##########################
#plt.figure()
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC1M1_0 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC1M1_1 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC1M2_0 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC1M2_1 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC1M3_0 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC1M3_1 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC2M1_0 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC2M1_1 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC2M2_0 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC2M2_1 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC2M3_0 for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.uvC2M3_1 for p in protos] )
#plt.title('cameras meas')
#
#
#plt.figure()
#plt.plot( [p.measurementsLastLatest.measurementsX.wIMU_0 for p in protos] )
#plt.plot( [p.measurementsLastLatest.measurementsX.wIMU_1 for p in protos] )
#plt.plot( [p.measurementsLastLatest.measurementsX.wIMU_2 for p in protos] )
#plt.plot( [p.measurementsLastLatest.measurementsX.aIMU_0 for p in protos] )
#plt.plot( [p.measurementsLastLatest.measurementsX.aIMU_1 for p in protos] )
#plt.plot( [p.measurementsLastLatest.measurementsX.aIMU_2 for p in protos] )
#plt.legend(['wx','wy','wz','ax','ay','az'])
#plt.title('imu meas')
#
#
#plt.figure()
#plt.plot( [p.measurementsCurrentX.measurementsX.ur for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsX.up for p in protos] )
#plt.legend(['ur','up'])
#plt.title('actuator meas')
#
#plt.figure()
#plt.plot( [p.measurementsCurrentX.measurementsU.dddelta for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsU.dur for p in protos] )
#plt.plot( [p.measurementsCurrentX.measurementsU.dup for p in protos] )
#plt.legend(['dddelta','dur','dup'])
#plt.title('dactuator meas')

plt.show()
