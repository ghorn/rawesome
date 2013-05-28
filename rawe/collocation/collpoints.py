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

import scipy.special
import numpy as np

def mkCollocationPoints(collPoly,deg):
#    assert(deg>0)
    # Legendre collocation points
    legendre_points = [[0],
                       [0,0.500000],
                       [0,0.211325,0.788675],
                       [0,0.112702,0.500000,0.887298],
                       [0,0.069432,0.330009,0.669991,0.930568],
                       [0,0.046910,0.230765,0.500000,0.769235,0.953090]]
    
    # Radau collocation points
    radau_points = [[0],
                    [0,1.000000],
                    [0,0.333333,1.000000],
                    [0,0.155051,0.644949,1.000000],
                    [0,0.088588,0.409467,0.787659,1.000000],
                    [0,0.057104,0.276843,0.583590,0.860240,1.000000]]

    hardCoded = None
    if collPoly == 'LEGENDRE':
        alpha = 0.0
        beta = 0.0
        roots = 0.5*(scipy.special.j_roots(deg,alpha,beta)[0] + 1)
        if deg<6:
            hardCoded = legendre_points[deg]
    elif collPoly == 'RADAU':
        alpha = 1.0
        beta = 0.0
        if deg > 1:
            roots = np.append(0.5*(scipy.special.j_roots(deg-1,alpha,beta)[0] + 1),1.0)
        elif deg == 1:
            roots = np.array([1.0])
        else:
            roots = np.array([])
        if deg<6:
            hardCoded = radau_points[deg]
    else:
        raise ValueError('poly must be "LEGENDRE" or "RADAU", you tried to use: "'+str(collPoly)+"\"")

    # if value is hard-coded, make sure it matches
    if hardCoded is not None:
        err = list(np.array(hardCoded[1:]) - roots)
        for e in err:
            if np.abs(e) > 1e-6:
                msg = "roots of Gauss-Jacobi polynomial don't match old hard-coded values"
                msg += "\nhard-coded:            "+str(hardCoded[1:])
                msg += "\nscipy.special.j_roots: "+str(list(roots))
                raise ValueError(msg)

    def toReal(val):
        if np.imag(val) == 0:
            return np.real(val).item()
        else:
            raise Exception("collocation point has non-zero imaginary part")
    return [0.0]+[toReal(r) for r in list(roots)]

if __name__ == '__main__':
    print "LEGENDRE:"
    for k in range(1,10):
        print mkCollocationPoints('LEGENDRE',k)
    print "RADAU:"
    for k in range(1,10):
        print mkCollocationPoints('RADAU',k)
