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
import numpy

def setupModel(dae, conf):
    '''
    take the dae that has x/z/u/p added to it already and return
    the states added to it and return mass matrix and rhs of the dae residual
    '''
    m =  conf['mass']
    g = conf['g']

    j1 =  conf['j1']
    j31 = conf['j31']
    j2 =  conf['j2']
    j3 =  conf['j3']

    jCarousel = conf['jCarousel']
    cfric = conf['cfric']

    zt = conf['zt']
    rA = conf['rArm']

    e11 = dae['e11']
    e12 = dae['e12']
    e13 = dae['e13']

    e21 = dae['e21']
    e22 = dae['e22']
    e23 = dae['e23']

    e31 = dae['e31']
    e32 = dae['e32']
    e33 = dae['e33']

    x =   dae['x']
    y =   dae['y']
    z =   dae['z']

    dx  =  dae['dx']
    dy  =  dae['dy']
    dz  =  dae['dz']

    w1 =  dae['w_bn_b_x']
    w2 =  dae['w_bn_b_y']
    w3 =  dae['w_bn_b_z']

    ddelta = dae['ddelta']

    r = dae['r']
    dr = dae['dr']

    ddr = dae['ddr']

    tc = dae['motor_torque'] #Carousel motor torque

    # wind model
    def getWind():
        if 'wind_model' not in conf:
            return 0
        if conf['wind_model']['name'] == 'wind_shear':
            # use a logarithmic wind shear model
            # wind(z) = w0 * log((z+zt)/zt) / log(z0/zt)
            # where w0 is wind at z0 altitude
            # zt is surface roughness characteristic length
            z0 = conf['wind_model']['z0']
            zt_roughness = conf['wind_model']['zt_roughness']
            zsat = 0.5*(-z+C.sqrt(z*z))
            return dae['w0']*C.log((zsat+zt_roughness+2)/zt_roughness)/C.log(z0/zt_roughness)
        elif conf['wind_model']['name'] == 'constant':
            # constant wind
            return dae['w0']
    wind_x = getWind()
    dae['wind_at_altitude'] = wind_x

    dp_carousel_frame = C.veccat( [ dx - ddelta*y
                                  , dy + ddelta*(rA + x)
                                  , dz
                                  ]) - C.veccat([dae['cos_delta']*wind_x, dae['sin_delta']*wind_x, 0])
    R_c2b = C.veccat( [dae[n] for n in ['e11', 'e12', 'e13',
                                        'e21', 'e22', 'e23',
                                        'e31', 'e32', 'e33']]
                      ).reshape((3,3))

    # Aircraft velocity w.r.t. inertial frame, given in its own reference frame
    # (needed to compute the aero forces and torques !)
    dpE = C.mul( R_c2b, dp_carousel_frame )

    (f1, f2, f3, t1, t2, t3) = aeroForcesTorques(dae, conf, dp_carousel_frame, dpE,
                                                 dae['w_bn_b'],
                                                 (dae['e21'], dae['e22'], dae['e23'])
                                                 )

    # mass matrix
    mm = C.SXMatrix(8, 8)
    mm[0,0] = jCarousel + m*rA*rA + m*x*x + m*y*y + 2*m*rA*x 
    mm[0,1] = -m*y 
    mm[0,2] = m*(rA + x) 
    mm[0,3] = 0 
    mm[0,4] = 0 
    mm[0,5] = 0 
    mm[0,6] = 0 
    mm[0,7] = 0
    
    mm[1,0] = -m*y 
    mm[1,1] = m 
    mm[1,2] = 0 
    mm[1,3] = 0 
    mm[1,4] = 0 
    mm[1,5] = 0 
    mm[1,6] = 0 
    mm[1,7] = x + zt*e31
    
    mm[2,0] = m*(rA + x) 
    mm[2,1] = 0 
    mm[2,2] = m 
    mm[2,3] = 0 
    mm[2,4] = 0 
    mm[2,5] = 0 
    mm[2,6] = 0 
    mm[2,7] = y + zt*e32
    
    mm[3,0] = 0 
    mm[3,1] = 0 
    mm[3,2] = 0 
    mm[3,3] = m 
    mm[3,4] = 0 
    mm[3,5] = 0 
    mm[3,6] = 0 
    mm[3,7] = z + zt*e33
    
    mm[4,0] = 0 
    mm[4,1] = 0 
    mm[4,2] = 0 
    mm[4,3] = 0 
    mm[4,4] = j1 
    mm[4,5] = 0 
    mm[4,6] = j31 
    mm[4,7] = -zt*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33)
    
    mm[5,0] = 0 
    mm[5,1] = 0 
    mm[5,2] = 0 
    mm[5,3] = 0 
    mm[5,4] = 0 
    mm[5,5] = j2 
    mm[5,6] = 0 
    mm[5,7] = zt*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33)
    
    mm[6,0] = 0 
    mm[6,1] = 0 
    mm[6,2] = 0 
    mm[6,3] = 0 
    mm[6,4] = j31 
    mm[6,5] = 0 
    mm[6,6] = j3 
    mm[6,7] = 0
    
    mm[7,0] = -zt*(e11*e23*x - e13*e21*x + e12*e23*y - e13*e22*y + zt*e11*e23*e31 - zt*e13*e21*e31 + zt*e12*e23*e32 - zt*e13*e22*e32) 
    mm[7,1] = x + zt*e31 
    mm[7,2] = y + zt*e32 
    mm[7,3] = z + zt*e33 
    mm[7,4] = -zt*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33) 
    mm[7,5] = zt*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) 
    mm[7,6] = 0 
    mm[7,7] = 0

    # right hand side
    zt2 = zt*zt
    rhs = C.veccat(
          [ tc - cfric*ddelta - f1*y + f2*(rA + x) + dy*m*(dx - 2*ddelta*y) - dx*m*(dy + 2*ddelta*rA + 2*ddelta*x) 
          , f1 + ddelta*m*(dy + ddelta*rA + ddelta*x) + ddelta*dy*m 
          , f2 - ddelta*m*(dx - ddelta*y) - ddelta*dx*m 
          , f3 + g*m 
          , t1 - w2*(j3*w3 + j31*w1) + j2*w2*w3 
          , t2 + w1*(j3*w3 + j31*w1) - w3*(j1*w1 + j31*w3) 
          , t3 + w2*(j1*w1 + j31*w3) - j2*w1*w2
          , ddr*r-(zt*w1*(e11*x+e12*y+e13*z+zt*e11*e31+zt*e12*e32+zt*e13*e33)+zt*w2*(e21*x+e22*y+e23*z+zt*e21*e31+zt*e22*e32+zt*e23*e33))*(w3-ddelta*e33)-dx*(dx-zt*e21*(w1-ddelta*e13)+zt*e11*(w2-ddelta*e23))-dy*(dy-zt*e22*(w1-ddelta*e13)+zt*e12*(w2-ddelta*e23))-dz*(dz-zt*e23*(w1-ddelta*e13)+zt*e13*(w2-ddelta*e23))+dr*dr+(w1-ddelta*e13)*(e21*(zt*dx-zt2*e21*(w1-ddelta*e13)+zt2*e11*(w2-ddelta*e23))+e22*(zt*dy-zt2*e22*(w1-ddelta*e13)+zt2*e12*(w2-ddelta*e23))+zt*e23*(dz+zt*e13*w2-zt*e23*w1)+zt*e33*(w1*z+zt*e33*w1+ddelta*e11*x+ddelta*e12*y+zt*ddelta*e11*e31+zt*ddelta*e12*e32)+zt*e31*(x+zt*e31)*(w1-ddelta*e13)+zt*e32*(y+zt*e32)*(w1-ddelta*e13))-(w2-ddelta*e23)*(e11*(zt*dx-zt2*e21*(w1-ddelta*e13)+zt2*e11*(w2-ddelta*e23))+e12*(zt*dy-zt2*e22*(w1-ddelta*e13)+zt2*e12*(w2-ddelta*e23))+zt*e13*(dz+zt*e13*w2-zt*e23*w1)-zt*e33*(w2*z+zt*e33*w2+ddelta*e21*x+ddelta*e22*y+zt*ddelta*e21*e31+zt*ddelta*e22*e32)-zt*e31*(x+zt*e31)*(w2-ddelta*e23)-zt*e32*(y+zt*e32)*(w2-ddelta*e23))
          ] )
 
    dRexp = C.SXMatrix(3,3)

    dRexp[0,0] = e21*(w3 - ddelta*e33) - e31*(w2 - ddelta*e23) 
    dRexp[0,1] = e31*(w1 - ddelta*e13) - e11*(w3 - ddelta*e33) 
    dRexp[0,2] = e11*(w2 - ddelta*e23) - e21*(w1 - ddelta*e13) 

    dRexp[1,0] = e22*(w3 - ddelta*e33) - e32*(w2 - ddelta*e23) 
    dRexp[1,1] = e32*(w1 - ddelta*e13) - e12*(w3 - ddelta*e33) 
    dRexp[1,2] = e12*(w2 - ddelta*e23) - e22*(w1 - ddelta*e13) 

    dRexp[2,0] = e23*w3 - e33*w2 
    dRexp[2,1] = e33*w1 - e13*w3 
    dRexp[2,2] = e13*w2 - e23*w1
    
    c =(x + zt*e31)**2/2 + (y + zt*e32)**2/2 + (z + zt*e33)**2/2 - r**2/2
    
    cdot =dx*(x + zt*e31) + dy*(y + zt*e32) + dz*(z + zt*e33) + zt*(w2 - ddelta*e23)*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) - zt*(w1 - ddelta*e13)*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33) - r*dr

#    ddx = dae['ddx']
#    ddy = dae['ddy']
#    ddz = dae['ddz']
#    dw1 = dae['dw1']
#    dw2 = dae['dw2']
#    dddelta = dae['dddelta']
    ddx = dae.ddt('dx')
    ddy = dae.ddt('dy')
    ddz = dae.ddt('dz')
    dw1 = dae.ddt('w_bn_b_x')
    dw2 = dae.ddt('w_bn_b_y')
    dddelta = dae.ddt('ddelta')
    
    cddot = -(w1-ddelta*e13)*(zt*e23*(dz+zt*e13*w2-zt*e23*w1)+zt*e33*(w1*z+zt*e33*w1+ddelta*e11*x+ddelta*e12*y+zt*ddelta*e11*e31+zt*ddelta*e12*e32)+zt*e21*(dx+zt*e11*w2-zt*e21*w1-zt*ddelta*e11*e23+zt*ddelta*e13*e21)+zt*e22*(dy+zt*e12*w2-zt*e22*w1-zt*ddelta*e12*e23+zt*ddelta*e13*e22)+zt*e31*(x+zt*e31)*(w1-ddelta*e13)+zt*e32*(y+zt*e32)*(w1-ddelta*e13))+(w2-ddelta*e23)*(zt*e13*(dz+zt*e13*w2-zt*e23*w1)-zt*e33*(w2*z+zt*e33*w2+ddelta*e21*x+ddelta*e22*y+zt*ddelta*e21*e31+zt*ddelta*e22*e32)+zt*e11*(dx+zt*e11*w2-zt*e21*w1-zt*ddelta*e11*e23+zt*ddelta*e13*e21)+zt*e12*(dy+zt*e12*w2-zt*e22*w1-zt*ddelta*e12*e23+zt*ddelta*e13*e22)-zt*e31*(x+zt*e31)*(w2-ddelta*e23)-zt*e32*(y+zt*e32)*(w2-ddelta*e23))-ddr*r+(zt*w1*(e11*x+e12*y+e13*z+zt*e11*e31+zt*e12*e32+zt*e13*e33)+zt*w2*(e21*x+e22*y+e23*z+zt*e21*e31+zt*e22*e32+zt*e23*e33))*(w3-ddelta*e33)+dx*(dx+zt*e11*w2-zt*e21*w1-zt*ddelta*e11*e23+zt*ddelta*e13*e21)+dy*(dy+zt*e12*w2-zt*e22*w1-zt*ddelta*e12*e23+zt*ddelta*e13*e22)+dz*(dz+zt*e13*w2-zt*e23*w1)+ddx*(x+zt*e31)+ddy*(y+zt*e32)+ddz*(z+zt*e33)-dr*dr+zt*(dw2-dddelta*e23)*(e11*x+e12*y+e13*z+zt*e11*e31+zt*e12*e32+zt*e13*e33)-zt*(dw1-dddelta*e13)*(e21*x+e22*y+e23*z+zt*e21*e31+zt*e22*e32+zt*e23*e33)-zt*dddelta*(e11*e23*x-e13*e21*x+e12*e23*y-e13*e22*y+zt*e11*e23*e31-zt*e13*e21*e31+zt*e12*e23*e32-zt*e13*e22*e32)

#    cddot = (zt*w1*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) + zt*w2*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33))*(w3 - ddelta*e33) + dx*(dx + zt*e11*w2 - zt*e21*w1 - zt*ddelta*e11*e23 + zt*ddelta*e13*e21) + dy*(dy + zt*e12*w2 - zt*e22*w1 - zt*ddelta*e12*e23 + zt*ddelta*e13*e22) + dz*(dz + zt*e13*w2 - zt*e23*w1) + ddx*(x + zt*e31) + ddy*(y + zt*e32) + ddz*(z + zt*e33) - (w1 - ddelta*e13)*(e21*(zt*dx - zt**2*e21*(w1 - ddelta*e13) + zt**2*e11*(w2 - ddelta*e23)) + e22*(zt*dy - zt**2*e22*(w1 - ddelta*e13) + zt**2*e12*(w2 - ddelta*e23)) + zt*e33*(z*w1 + ddelta*e11*x + ddelta*e12*y + zt*e33*w1 + zt*ddelta*e11*e31 + zt*ddelta*e12*e32) + zt*e23*(dz + zt*e13*w2 - zt*e23*w1) + zt*e31*(w1 - ddelta*e13)*(x + zt*e31) + zt*e32*(w1 - ddelta*e13)*(y + zt*e32)) + (w2 - ddelta*e23)*(e11*(zt*dx - zt**2*e21*(w1 - ddelta*e13) + zt**2*e11*(w2 - ddelta*e23)) + e12*(zt*dy - zt**2*e22*(w1 - ddelta*e13) + zt**2*e12*(w2 - ddelta*e23)) - zt*e33*(z*w2 + ddelta*e21*x + ddelta*e22*y + zt*e33*w2 + zt*ddelta*e21*e31 + zt*ddelta*e22*e32) + zt*e13*(dz + zt*e13*w2 - zt*e23*w1) - zt*e31*(w2 - ddelta*e23)*(x + zt*e31) - zt*e32*(w2 - ddelta*e23)*(y + zt*e32)) + zt*(dw2 - dddelta*e23)*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) - zt*(dw1 - dddelta*e13)*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33) - zt*dddelta*(e11*e23*x - e13*e21*x + e12*e23*y - e13*e22*y + zt*e11*e23*e31 - zt*e13*e21*e31 + zt*e12*e23*e32 - zt*e13*e22*e32)
#      where
#        dw1 = dw @> 0
#        dw2 = dw @> 1
#        {-
#        dw3 = dw @> 2
#        -}
#        ddx = ddX @> 0
#        ddy = ddX @> 1
#        ddz = ddX @> 2
#        dddelta = dddelta' @> 0

    dae['c'] =  c
    dae['cdot'] = cdot
    dae['cddot'] = cddot
    return (mm, rhs, dRexp)
        
def carouselModel(conf):
    '''
    pass this a conf, and it'll return you a dae
    '''
    # empty Dae
    dae = Dae()
        
    # add some differential states/algebraic vars/controls/params
    dae.addZ("nu")
    dae.addX( [ "x"
              , "y"
              , "z"
              , "e11"
              , "e12"
              , "e13"
              , "e21"
              , "e22"
              , "e23"
              , "e31"
              , "e32"
              , "e33"
              , "dx"
              , "dy"
              , "dz"
              , "w_bn_b_x"
              , "w_bn_b_y"
              , "w_bn_b_z"
              , "ddelta"
              , "r"
              , "dr"
              , "aileron"
              , "elevator"
              , "motor_torque"
              , "ddr"
              ] )
    if conf['delta_parameterization'] == 'linear':
        dae.addX('delta')
        dae['cos_delta'] = C.cos(dae['delta'])
        dae['sin_delta'] = C.sin(dae['delta'])
        dae_delta_residual = dae.ddt('delta') - dae['ddelta']

    elif conf['delta_parameterization'] == 'cos_sin':
        dae.addX("cos_delta")
        dae.addX("sin_delta")
        norm = dae['cos_delta']**2 + dae['sin_delta']**2
        
        if 'stabilize_invariants' in conf and conf['stabilize_invariants'] == True:
            pole_delta = 0.5
        else:
            pole_delta = 0.0
        
        cos_delta_dot_st = -pole_delta/2.* ( dae['cos_delta'] - dae['cos_delta'] / norm )
        sin_delta_dot_st = -pole_delta/2.* ( dae['sin_delta'] - dae['sin_delta'] / norm )
        dae_delta_residual = C.veccat([dae.ddt('cos_delta') - (-dae['sin_delta']*dae['ddelta'] + cos_delta_dot_st),
                                       dae.ddt('sin_delta') - ( dae['cos_delta']*dae['ddelta'] + sin_delta_dot_st) ])
    else:
        raise ValueError('unrecognized delta_parameterization "'+conf['delta_parameterization']+'"')

    dae.addU( [ "daileron"
              , "delevator"
              , "dmotor_torque"
              , 'dddr'
              ] )
    # add wind parameter if wind shear is in configuration
    if 'wind_model' in conf:
        dae.addP( ['w0'] )

    # set some state derivatives as outputs
    dae['ddx'] = dae.ddt('dx')
    dae['ddy'] = dae.ddt('dy')
    dae['ddz'] = dae.ddt('dz')
    dae['ddt_w_bn_b_x'] = dae.ddt('w_bn_b_x')
    dae['ddt_w_bn_b_y'] = dae.ddt('w_bn_b_y')
    dae['ddt_w_bn_b_z'] = dae.ddt('w_bn_b_z')
    dae['w_bn_b'] = C.veccat([dae['w_bn_b_x'], dae['w_bn_b_y'], dae['w_bn_b_z']])

    # some outputs in for plotting
    dae['RPM'] = dae['ddelta']*60/(2*C.pi)
    dae['aileron_deg'] = dae['aileron']*180/C.pi
    dae['elevator_deg'] = dae['elevator']*180/C.pi
    dae['daileron_deg_s'] = dae['daileron']*180/C.pi
    dae['delevator_deg_s'] = dae['delevator']*180/C.pi
    
    dae['motor_power'] = dae['motor_torque']*dae['ddelta']

    dae['tether_tension'] = dae['r']*dae['nu']
    dae['winch_power'] = -dae['tether_tension']*dae['dr']
    
    dae['dcm'] = C.vertcat([C.horzcat([dae['e11'],dae['e12'],dae['e13']]),
                            C.horzcat([dae['e21'],dae['e22'],dae['e23']]),
                            C.horzcat([dae['e31'],dae['e32'],dae['e33']])])
    
    # line angle
    dae['cos_line_angle'] = \
      -(dae['e31']*dae['x'] + dae['e32']*dae['y'] + dae['e33']*dae['z']) / C.sqrt(dae['x']**2 + dae['y']**2 + dae['z']**2)
    dae['line_angle_deg'] = C.arccos(dae['cos_line_angle'])*180.0/C.pi

    (massMatrix, rhs, dRexp) = setupModel(dae, conf)
    
    if 'stabilize_invariants' in conf and conf['stabilize_invariants'] == True:
        RotPole = 0.5
    else:
        RotPole = 0.0
    Rst = RotPole*C.mul( dae['dcm'], (C.inv(C.mul(dae['dcm'].T,dae['dcm'])) - numpy.eye(3)) )
    
    ode = C.veccat([
        C.veccat([dae.ddt(name) for name in ['x','y','z']]) - C.veccat([dae['dx'],dae['dy'],dae['dz']]),
        C.veccat([dae.ddt(name) for name in ["e11","e12","e13",
                                             "e21","e22","e23",
                                             "e31","e32","e33"]]) - ( dRexp.trans().reshape([9,1]) + Rst.reshape([9,1]) ),
        dae_delta_residual,
#        C.veccat([dae.ddt(name) for name in ['dx','dy','dz']]) - C.veccat([dae['ddx'],dae['ddy'],dae['ddz']]),
#        C.veccat([dae.ddt(name) for name in ['w1','w2','w3']]) - C.veccat([dae['dw1'],dae['dw2'],dae['dw3']]),
#        dae.ddt('ddelta') - dae['dddelta'],
        dae.ddt('r') - dae['dr'],
        dae.ddt('dr') - dae['ddr'],
        dae.ddt('aileron') - dae['daileron'],
        dae.ddt('elevator') - dae['delevator'],
        dae.ddt('motor_torque') - dae['dmotor_torque'],
        dae.ddt('ddr') - dae['dddr']
        ])

    if 'stabilize_invariants' in conf and conf['stabilize_invariants'] == True:
        cPole = 0.5
    else:
        cPole = 0.0
    rhs[-1] -= 2*cPole*dae['cdot'] + cPole*cPole*dae['c']
    
    psuedoZVec = C.veccat([dae.ddt(name) for name in ['ddelta','dx','dy','dz','w_bn_b_x','w_bn_b_y','w_bn_b_z']]+[dae['nu']])
    alg = C.mul(massMatrix, psuedoZVec) - rhs
    dae.setResidual([ode,alg])
    
    return dae
