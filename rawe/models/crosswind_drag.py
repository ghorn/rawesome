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
from rawe.dae import Dae
from aero import aeroForcesTorques

def get_wind(dae, conf):
    '''
    return the wind at current altitude based on configuration
    '''
    if conf['wind_model']['name'] == 'wind_shear':
        # use a logarithmic wind shear model
        # wind(z) = w0 * log((altitude+zt)/zt) / log(z0/zt)
        # where w0 is wind at z0 altitude
        # zt is surface roughness characteristic length
        # altitude is -z - altitude0, where altitude0 is a user parameter
        z = dae['r_n2b_n_z']
        z0 = conf['wind_model']['z0']
        zt_roughness = conf['wind_model']['zt_roughness']
        altitude = -z - conf['wind_model']['altitude0']
        return dae['w0']*C.log((altitude+zt_roughness)/zt_roughness) / \
                         C.log(z0/zt_roughness)
    elif conf['wind_model']['name'] == 'constant':
        # constant wind
        return dae['w0']

def skew_mat(xyz):
    '''
    return the skew symmetrix matrix of a 3-vector
    '''
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    return C.vertcat([C.horzcat([ 0, -z,  y]),
                      C.horzcat([ z,  0, -x]),
                      C.horzcat([-y,  x,  0])])

def compute_mass_matrix(dae, conf, f1, f2, f3, t1, t2, t3):
    '''
    take the dae that has x/z/u/p added to it already and return
    the states added to it and return mass matrix and rhs of the dae residual
    '''

    R_b2n = dae['R_n2b'].T
    r_n2b_n = C.veccat([dae['r_n2b_n_x'], dae['r_n2b_n_y'], dae['r_n2b_n_z']])
    r_b2bridle_b = C.veccat([0,0,conf['zt']])
    v_bn_n = C.veccat([dae['v_bn_n_x'],dae['v_bn_n_y'],dae['v_bn_n_z']])

    r_n2bridle_n = r_n2b_n + C.mul(R_b2n, r_b2bridle_b)

    mm00 = C.diag([1,1,1]) * (conf['mass'] + conf['tether_mass']/3.0)
    mm01 = C.SXMatrix(3,3)
    mm10 = mm01.T
    mm02 = r_n2bridle_n
    mm20 = mm02.T
    J = C.vertcat([C.horzcat([ conf['j1'],          0, conf['j31']]),
                   C.horzcat([          0, conf['j2'],           0]),
                   C.horzcat([conf['j31'],          0,  conf['j3']])])
    mm11 = J
    mm12 = C.cross(r_b2bridle_b, C.mul(dae['R_n2b'], r_n2b_n))
    mm21 = mm12.T
    mm22 = C.SXMatrix(1,1)

    mm = C.vertcat([C.horzcat([mm00,mm01,mm02]),
                    C.horzcat([mm10,mm11,mm12]),
                    C.horzcat([mm20,mm21,mm22])])

    # right hand side
    rhs0 = C.veccat([f1,f2,f3 + conf['g']*(conf['mass'] + conf['tether_mass']*0.5)])
    rhs1 = C.veccat([t1,t2,t3]) - C.cross(dae['w_bn_b'], C.mul(J, dae['w_bn_b']))

    # last element of RHS
    R_n2b = dae['R_n2b']
    w_bn_b = dae['w_bn_b']
    grad_r_cdot = v_bn_n + C.mul(R_b2n, C.cross(dae['w_bn_b'], r_b2bridle_b))
    tPR = - C.cross(C.mul(R_n2b, r_n2b_n), C.cross(w_bn_b, r_b2bridle_b)) - \
          C.cross(C.mul(R_n2b, v_bn_n), r_b2bridle_b)
    rhs2 = -C.mul(grad_r_cdot.T, v_bn_n) - C.mul(tPR.T, w_bn_b)

    rhs = C.veccat([rhs0,rhs1,rhs2])

    c = 0.5*(C.mul(r_n2bridle_n.T, r_n2bridle_n) - dae['r']**2)
    v_bridlen_n = v_bn_n + C.mul(R_b2n, C.cross(w_bn_b, r_b2bridle_b))
    cdot = C.mul(r_n2bridle_n.T, v_bridlen_n)

    a_bn_n = C.veccat([dae.ddt(name) for name in ['v_bn_n_x','v_bn_n_y','v_bn_n_z']])
    dw_bn_b = C.veccat([dae.ddt(name) for name in ['w_bn_b_x','w_bn_b_y','w_bn_b_z']])
    a_bridlen_n = a_bn_n + C.mul(R_b2n, C.cross(dw_bn_b, r_b2bridle_b) + C.cross(w_bn_b, C.cross(w_bn_b, r_b2bridle_b)))
    cddot = C.mul(v_bridlen_n.T, v_bridlen_n) + C.mul(r_n2bridle_n.T, a_bridlen_n)

    dae['c'] = c
    dae['cdot'] = cdot
    dae['cddot'] = cddot

    return (mm, rhs)

def crosswind_model(conf):
    '''
    pass this a conf, and it'll return you a dae
    '''
    # empty Dae
    dae = Dae()

    # add some differential states/algebraic vars/controls/params
    dae.addZ('nu')
    dae.addX( ['r_n2b_n_x', 'r_n2b_n_y', 'r_n2b_n_z',
               'e11', 'e12', 'e13',
               'e21', 'e22', 'e23',
               'e31', 'e32', 'e33',
               'v_bn_n_x', 'v_bn_n_y', 'v_bn_n_z',
               'w_bn_b_x', 'w_bn_b_y', 'w_bn_b_z',
               'aileron', 'elevator', 'rudder', 'flaps', 'prop_drag'
           ])
    dae.addU( [ 'daileron'
              , 'delevator'
              , 'drudder'
              , 'dflaps'
              , 'dprop_drag'
              ] )
    dae.addP( ['r','w0'] )

    # set some state derivatives as outputs
    dae['ddx'] = dae.ddt('v_bn_n_x')
    dae['ddy'] = dae.ddt('v_bn_n_y')
    dae['ddz'] = dae.ddt('v_bn_n_z')
    dae['ddt_w_bn_b_x'] = dae.ddt('w_bn_b_x')
    dae['ddt_w_bn_b_y'] = dae.ddt('w_bn_b_y')
    dae['ddt_w_bn_b_z'] = dae.ddt('w_bn_b_z')
    dae['w_bn_b'] = C.veccat([dae['w_bn_b_x'], dae['w_bn_b_y'], dae['w_bn_b_z']])

    # some outputs in degrees for plotting
    dae['aileron_deg']     = dae['aileron']*180/C.pi
    dae['elevator_deg']    = dae['elevator']*180/C.pi
    dae['rudder_deg']      = dae['rudder']*180/C.pi
    dae['daileron_deg_s']  = dae['daileron']*180/C.pi
    dae['delevator_deg_s'] = dae['delevator']*180/C.pi
    dae['drudder_deg_s']   = dae['drudder']*180/C.pi

    # tether tension == radius*nu where nu is alg. var associated with x^2+y^2+z^2-r^2==0
    dae['tether_tension'] = dae['r']*dae['nu']

    dae['R_n2b'] = C.vertcat([C.horzcat([dae['e11'], dae['e12'], dae['e13']]),
                              C.horzcat([dae['e21'], dae['e22'], dae['e23']]),
                              C.horzcat([dae['e31'], dae['e32'], dae['e33']])])

    # local wind
    dae['wind_at_altitude'] = get_wind(dae, conf)

    # wind velocity in NED
    dae['v_wn_n'] = C.veccat([dae['wind_at_altitude'], 0, 0])

    # body velocity in NED
    dae['v_bn_n'] = C.veccat( [ dae['v_bn_n_x'] , dae['v_bn_n_y'], dae['v_bn_n_z'] ] )

    # body velocity w.r.t. wind in NED
    v_bw_n =  dae['v_bn_n'] - dae['v_wn_n']

    # body velocity w.r.t. wind in body frame
    v_bw_b = C.mul( dae['R_n2b'], v_bw_n )

    # compute aerodynamic forces and moments
    (f1, f2, f3, t1, t2, t3) = \
        aeroForcesTorques(dae, conf, v_bw_n, v_bw_b,
                          dae['w_bn_b'],
                          (dae['e21'], dae['e22'], dae['e23']))
    dae['prop_drag_vector'] = dae['prop_drag']*v_bw_n/dae['airspeed']
    dae['prop_power'] = C.mul(dae['prop_drag_vector'].T, v_bw_n)
    f1 -= dae['prop_drag_vector'][0]
    f2 -= dae['prop_drag_vector'][1]
    f3 -= dae['prop_drag_vector'][2]

    # if we are running a homotopy,
    # add psudeo forces and moments as algebraic states
    if 'runHomotopy' in conf and conf['runHomotopy']:
        gamma_homotopy = dae.addP('gamma_homotopy')
        f1 = f1*gamma_homotopy + dae.addZ('f1_homotopy')*(1 - gamma_homotopy)
        f2 = f2*gamma_homotopy + dae.addZ('f2_homotopy')*(1 - gamma_homotopy)
        f3 = f3*gamma_homotopy + dae.addZ('f3_homotopy')*(1 - gamma_homotopy)
        t1 = t1*gamma_homotopy + dae.addZ('t1_homotopy')*(1 - gamma_homotopy)
        t2 = t2*gamma_homotopy + dae.addZ('t2_homotopy')*(1 - gamma_homotopy)
        t3 = t3*gamma_homotopy + dae.addZ('t3_homotopy')*(1 - gamma_homotopy)

    # derivative of dcm
    dRexp = C.mul(skew_mat(dae['w_bn_b']).T, dae['R_n2b'])
    ddt_R_n2b = C.vertcat(\
        [C.horzcat([dae.ddt(name) for name in ['e11', 'e12', 'e13']]),
         C.horzcat([dae.ddt(name) for name in ['e21', 'e22', 'e23']]),
         C.horzcat([dae.ddt(name) for name in ['e31', 'e32', 'e33']])])

    # get mass matrix, rhs
    (mass_matrix, rhs) = compute_mass_matrix(dae, conf, f1, f2, f3, t1, t2, t3)

    # set up the residual
    ode = C.veccat([
        C.veccat([dae.ddt('r_n2b_n_'+name) - dae['v_bn_n_'+name] for name in ['x','y','z']]),
        C.veccat([ddt_R_n2b - dRexp]),
        dae.ddt('aileron') - dae['daileron'],
        dae.ddt('elevator') - dae['delevator'],
        dae.ddt('rudder') - dae['drudder'],
        dae.ddt('flaps') - dae['dflaps'],
        dae.ddt('prop_drag') - dae['dprop_drag']
        ])

    # acceleration for plotting
    ddx = dae['ddx']
    ddy = dae['ddy']
    ddz = dae['ddz']
    dae['accel_g'] = C.sqrt(ddx**2 + ddy**2 + (ddz + 9.8)**2)/9.8
    dae['accel_without_gravity_g'] = C.sqrt(ddx**2 + ddy**2 + ddz**2)/9.8
    dae['accel'] = C.sqrt(ddx**2 + ddy**2 + (ddz+9.8)**2)
    dae['accel_without_gravity'] = C.sqrt(ddx**2 + ddy**2 + ddz**2)

    # line angle
    dae['cos_line_angle'] = \
      -(dae['e31']*dae['r_n2b_n_x'] + dae['e32']*dae['r_n2b_n_y'] + dae['e33']*dae['r_n2b_n_z']) / \
       C.sqrt(dae['r_n2b_n_x']**2 + dae['r_n2b_n_y']**2 + dae['r_n2b_n_z']**2)
    dae['line_angle_deg'] = C.arccos(dae['cos_line_angle'])*180.0/C.pi

    # add local loyd's limit
    def addLoydsLimit():
        w = dae['wind_at_altitude']
        cL = dae['cL']
        cD = dae['cD'] + dae['cD_tether']
        rho = conf['rho']
        S = conf['sref']
        loyds = 2/27.0*rho*S*w**3*cL**3/cD**2
        loyds2 = 2/27.0*rho*S*w**3*(cL**2/cD**2 + 1)**(1.5)*cD
        dae["loyds_limit"] = loyds
        dae["loyds_limit_exact"] = loyds2
    addLoydsLimit()

    psuedo_zvec = C.veccat([dae.ddt(name) for name in \
        ['v_bn_n_x','v_bn_n_y','v_bn_n_z','w_bn_b_x','w_bn_b_y','w_bn_b_z']]+[dae['nu']])
    alg = C.mul(mass_matrix, psuedo_zvec) - rhs
    dae.setResidual( [ode, alg] )

    return dae
