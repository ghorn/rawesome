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

from casadi import pi

import rawe
from carouselSteadyState import getSteadyState

if __name__=='__main__':
    from conf import conf
    
    print "creating model..."
    dae = rawe.models.carousel(conf)
    dae.convertToOde()

    name = 'carouselOde'
    # write the file that computes f,MM in 0 == f(x,u,p) + MM(x,u,p)*[xdot;z]
    blah = dae.octaveSimGen(name)
    f = open(name+'_modelAndJacob.m','w')
    f.write(blah)
    f.close()

    # write the file that computes xdot,z from f,MM
    f = open(name+'_xDotAndZ.m','w')
    f.write('''\
function [xDot,z] = %(name)s_xDotAndZ(x,u,p)
[f,MM] = %(name)s_modelAndJacob(x,u,p);
xDotAndZ = -(MM\\f);
xDot = xDotAndZ(1:%(nx)d);
z = xDotAndZ(%(nx)d+1:end);
end
''' % {'name':name,'nx':len(dae.xNames())})
    f.close()

    steadyState = getSteadyState(dae,conf,2*pi,1.2)
    x = {}
    lines = []
    lines.append('function [x,z,u,p] = carouselOde_steadyState()')
    lines.append('x = zeros('+str(len(dae.xNames()))+',1);')
    lines.append('u = zeros('+str(len(dae.uNames()))+',1);')
    lines.append('p = zeros('+str(len(dae.pNames()))+',1);')
    for k,name in enumerate(dae.xNames()):
        lines.append('x('+str(k+1)+') = '+str(steadyState[name])+'; % '+name)
    for k,name in enumerate(dae.zNames()):
        lines.append('z('+str(k+1)+') = '+str(steadyState[name])+'; % '+name)
    for k,name in enumerate(dae.uNames()):
        lines.append('u('+str(k+1)+') = '+str(steadyState[name])+'; % '+name)
    for k,name in enumerate(dae.pNames()):
        lines.append('p('+str(k+1)+') = '+str(steadyState[name])+'; % '+name)
    lines.append('end')

    f = open('carouselOde_steadyState.m','w')
    f.write('\n'.join(lines))
    f.close()
