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

import rawe
from rawe.models import arianne_conf
import casadi as C
import numpy as np
#import os


def makeDae( conf = None ):
    if conf is None:
        conf = arianne_conf.makeConf()

    # Make model
    dae = rawe.models.carousel(conf)

    (xDotSol, zSol) = dae.solveForXDotAndZ()

    # Get variables and outputs from the model
    ddp = C.vertcat([xDotSol['dx'],xDotSol['dy'],xDotSol['dz']])
    ddt_w_bn_b  = C.vertcat([xDotSol['w_bn_b_x'],xDotSol['w_bn_b_y'],xDotSol['w_bn_b_z']])
    x =   dae['x']
    y =   dae['y']
    z =   dae['z']

    e31 = dae['e31']
    e32 = dae['e32']
    e33 = dae['e33']

    zt = conf['zt']

    dx  =  dae['dx']
    dy  =  dae['dy']

    ddelta = dae['ddelta']
    dddelta = xDotSol['ddelta']

    R = dae['R_c2b']
    rA = conf['rArm']
    g = conf['g']

    # Rotation matrix to convert from NWU to NED frame type
    R_nwu2ned = np.eye( 3 )
    R_nwu2ned[1, 1] = R_nwu2ned[2, 2] = -1.0

    ############################################################################
    #
    # IMU model
    #
    ############################################################################

    # Load IMU position and orientation w.r.t. body frame
#    pIMU = C.mul(R_nwu2ned, C.DMatrix(np.loadtxt(os.path.join(propertiesDir,'IMU/pIMU.dat'))))
    pIMU = C.DMatrix([0,0,0])
    pIMU = C.mul(R_nwu2ned, pIMU)
#0
#0
#0

#    RIMU = C.mul(R_nwu2ned, C.DMatrix(np.loadtxt(os.path.join(propertiesDir,'IMU/RIMU.dat'))))

#9.937680e-01   6.949103e-02    8.715574e-02
#-6.975647e-02  9.975641e-01    0
#-8.694344e-02  -6.079677e-03   9.961947e-01
    RIMU = C.DMatrix([[9.937680e-01,   6.949103e-02, 8.715574e-02],
                      [-6.975647e-02,  9.975641e-01,            0],
                      [-8.694344e-02, -6.079677e-03, 9.961947e-01]])
    RIMU = C.mul(R_nwu2ned, RIMU)

        # Define IMU measurement functions
        # TODO here is omitted the term: w x w pIMU
    # The sign of gravity is negative because of the NED convention (z points down!)
    ddpIMU_c = ddp - ddelta ** 2 * C.vertcat([x + rA, y, 0]) + 2 * ddelta * C.vertcat([-dy, dx, 0]) + \
                dddelta * C.vertcat([-y, x + rA, 0]) - C.vertcat([0, 0, g])
    ddpIMU = C.mul(R, ddpIMU_c)
    aBridle = C.cross(ddt_w_bn_b, pIMU)
    ddpIMU += aBridle
    ddpIMU = C.mul(RIMU,ddpIMU)
    # You can add a parameter to conf which will give the model 3 extra states with derivative 0 (to act as parameter) for the bias in the acceleration measurements. If that is present, it should be added to the measurement of the acceleration
    if 'useIMUAccelerationBias' in conf and conf['useIMUAccelerationBias']:
        IMUAccelerationBias = C.vertcat([dae['IMUAccelerationBias1'],dae['IMUAccelerationBias2'],dae['IMUAccelerationBias3']])
        ddpIMU += IMUAccelerationBias

    # For the accelerometers
    dae['IMU_acceleration'] = ddpIMU
    dae['IMU_acceleration_x'] = dae['IMU_acceleration'][0]
    dae['IMU_acceleration_y'] = dae['IMU_acceleration'][1]
    dae['IMU_acceleration_z'] = dae['IMU_acceleration'][2]

    # ... and for the gyroscopes
    dae['IMU_angular_velocity'] = C.mul(RIMU, dae['w_bn_b'])
    dae['IMU_angular_velocity_x'] = dae['IMU_angular_velocity'][0]
    dae['IMU_angular_velocity_y'] = dae['IMU_angular_velocity'][1]
    dae['IMU_angular_velocity_z'] = dae['IMU_angular_velocity'][2]

    if 'kinematicIMUAccelerationModel' in conf and conf['kinematicIMUAccelerationModel']:
        dae['vel_error_x'] = dae['dx']-dae['dx_IMU']
        dae['vel_error_y'] = dae['dy']-dae['dy_IMU']
        dae['vel_error_z'] = dae['dz']-dae['dz_IMU']

        ############################################################################
        #
        # LAS
        #
        ############################################################################

    x_tether = (x + zt*e31)
    y_tether = (y + zt*e32)
    z_tether = (z + zt*e33)

    dae['lineAngles'] = C.vertcat([C.arctan(y_tether/x_tether),C.arctan(z_tether/x_tether)])
    dae['lineAngle_hor'] = dae['lineAngles'][0]
    dae['lineAngle_ver'] = dae['lineAngles'][1]

    return dae
