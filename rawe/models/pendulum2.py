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

from rawe.dae import Dae
import casadi as C

def pendulumModel(nSteps=None):
    dae = Dae()
    dae.addZ( [ "tau"
              ] )
    dae.addX( [ "x"
              , "z"
              , "dx"
              , "dz"
              ] )
    dae.addU( [ "torque"
              ] )
    dae.addP( [ "m"
              ] )

    r = 0.3
    
    ddx = dae.ddt('dx')
    ddz = dae.ddt('dz')

    ode = [ dae['dx'] - dae.ddt('x')
          , dae['dz'] - dae.ddt('z')
#          , dae['ddx'] - dae.ddt('dx')
#          , dae['ddz'] - dae.ddt('dz')
          ]

    fx =  dae['torque']*dae['z']
    fz = -dae['torque']*dae['x'] + dae['m']*9.8
    
    alg = [ dae['m']*ddx + dae['x']*dae['tau'] - fx
          , dae['m']*ddz + dae['z']*dae['tau'] - fz
          , dae['x']*ddx + dae['z']*ddz + (dae['dx']*dae['dx'] + dae['dz']*dae['dz']) ]

    dae['c']    = dae['x']*dae['x'] + dae['z']*dae['z'] - r*r
    dae['cdot'] = dae['dx']*dae['x'] + dae['dz']*dae['z']

    dae.setResidual( [ode, alg] )

    return dae

if __name__=='__main__':
    pendulumModel()
    pendulumModel(nSteps=10)
