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

def makeDae():
    dae = rawe.Dae()

    [pos,vel,mass] = dae.addX( ["pos","vel","mass"] )
    thrust = dae.addU( "thrust" )
    
    # some extra outputs, why not
    dae['posvel'] = pos*vel
    dae['velvel'] = vel*vel

    # residual
    dae.setResidual([dae.ddt('pos') - vel,
                     dae.ddt('vel') - thrust/mass,
                     dae.ddt('mass') - 0.1*thrust*thrust])

    return dae
