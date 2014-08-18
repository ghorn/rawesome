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

def getEuler(ocp, k):
    r11 = ocp.lookup('e11',timestep=k)
    r12 = ocp.lookup('e12',timestep=k)
    mr13 = -ocp.lookup('e13',timestep=k)
#     mr13 -- nan protect
#       | mr13' >  1 =  1
#       | mr13' < -1 = -1
#       | otherwise = mr13'
    r23 = ocp.lookup('e23',timestep=k)
    r33 = ocp.lookup('e33',timestep=k)

    yaw   = C.arctan2(r12,r11)
    pitch = C.arcsin(mr13)
    roll  = C.arctan2(r23,r33)
    return (yaw,pitch,roll)

# euler angle periodic constraints
def periodicEulers(ocp):
    (yaw0,pitch0,roll0) = getEuler(ocp, 0)
    (yawF,pitchF,rollF) = getEuler(ocp, -1)
    ocp.constrain(yaw0,'==',yawF)
    ocp.constrain(pitch0,'==',pitchF)
    ocp.constrain(roll0,'==',rollF)


def getDcm(ocp,k,prefix='e'):
    m11 = ocp.lookup(prefix+'11',timestep=k)
    m12 = ocp.lookup(prefix+'12',timestep=k)
    m13 = ocp.lookup(prefix+'13',timestep=k)

    m21 = ocp.lookup(prefix+'21',timestep=k)
    m22 = ocp.lookup(prefix+'22',timestep=k)
    m23 = ocp.lookup(prefix+'23',timestep=k)

    m31 = ocp.lookup(prefix+'31',timestep=k)
    m32 = ocp.lookup(prefix+'32',timestep=k)
    m33 = ocp.lookup(prefix+'33',timestep=k)

    return C.vertcat([C.horzcat([m11,m12,m13]),
                      C.horzcat([m21,m22,m23]),
                      C.horzcat([m31,m32,m33])])

def getOrthonormalizedDcm(ocp,k):
    m = getDcm(ocp,k)
    return orthonormalizeDcm(m)

def get_orthonormal_constraints(R):
    ret = []
    ret.append( (C.mul(R[0,:],R[0,:].T) - 1,  'R1[0]: e1^T * e1 - 1 == 0') )
    ret.append( (C.mul(R[1,:],R[0,:].T),      'R1[0]: e2^T * e1     == 0') )
    ret.append( (C.mul(R[1,:],R[1,:].T) - 1,  'R1[0]: e2^T * e2 - 1 == 0') )
    rhon = C.cross(R[0,:],R[1,:]) - R[2,:]
    ret.append( (rhon[0,0], 'R1[0]: ( e1^T X e2 - e3 )[0] == 0') )
    ret.append( (rhon[0,2], 'R1[0]: ( e1^T X e2 - e3 )[1] == 0') )
    ret.append( (rhon[0,1], 'R1[0]: ( e1^T X e2 - e3 )[2] == 0') )
    return ret

def makeOrthonormal(ocp_,R):
    for val,description in get_orthonormal_constraints(R):
        ocp_.constrain(val, '==', 0, tag=(description, None))

def matchDcms(ocp,R0,Rf,tag=None):
    err = C.mul(R0.T, Rf)
    if tag is None:
        tag = ''
    else:
        tag = ' '+tag
    ocp.constrain(err[0,1], '==', 0, tag=('dcm matching 01'+tag,None))
    ocp.constrain(err[0,2], '==', 0, tag=('dcm matching 02'+tag,None))
    ocp.constrain(err[1,2], '==', 0, tag=('dcm matching 12'+tag,None))

    ocp.constrain(err[0,0], '>=', 0.5, tag=('dcm matching 00'+tag,None))
    ocp.constrain(err[1,1], '>=', 0.5, tag=('dcm matching 11'+tag,None))
    ocp.constrain(err[2,2], '>=', 0.5, tag=('dcm matching 22'+tag,None))

def periodicDcm(ocp):
    R0 = getDcm(ocp,0)
    Rf = getDcm(ocp,-1)
    matchDcms(ocp,R0,Rf)

# dcm periodic constraints
def periodicOrthonormalizedDcm(ocp):
    R0 = getOrthonormalizedDcm(ocp,0)
    Rf = getOrthonormalizedDcm(ocp,-1)
    matchDcms(ocp,R0,Rf)

def orthonormalizeDcm(m):
    ## OGRE (www.ogre3d.org) is made available under the MIT License.
    ##
    ## Copyright (c) 2000-2009 Torus Knot Software Ltd
    ##
    ## Permission is hereby granted, free of charge, to any person obtaining a copy
    ## of this software and associated documentation files (the "Software"), to deal
    ## in the Software without restriction, including without limitation the rights
    ## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    ## copies of the Software, and to permit persons to whom the Software is
    ## furnished to do so, subject to the following conditions:
    ##
    ## The above copyright notice and this permission notice shall be included in
    ## all copies or substantial portions of the Software.
    ##
    ## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    ## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    ## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    ## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    ## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    ## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    ## THE SOFTWARE.

    # Algorithm uses Gram-Schmidt orthogonalization.  If 'this' matrix is
    # M = [m0|m1|m2], then orthonormal output matrix is Q = [q0|q1|q2],
    #
    #   q0 = m0/|m0|
    #   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
    #   q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|
    #
    # where |V| indicates length of vector V and A*B indicates dot
    # product of vectors A and B.

    m00 = m[0,0]
    m01 = m[0,1]
    m02 = m[0,2]

    m10 = m[1,0]
    m11 = m[1,1]
    m12 = m[1,2]

    m20 = m[2,0]
    m21 = m[2,1]
    m22 = m[2,2]

    # compute q0
    fInvLength = 1.0/C.sqrt(m00*m00 + m10*m10 + m20*m20)

    m00 *= fInvLength
    m10 *= fInvLength
    m20 *= fInvLength

    # compute q1
    fDot0 = m00*m01 + m10*m11 + m20*m21

    m01 -= fDot0*m00
    m11 -= fDot0*m10
    m21 -= fDot0*m20

    fInvLength = 1.0/C.sqrt(m01*m01 + m11*m11 + m21*m21)

    m01 *= fInvLength
    m11 *= fInvLength
    m21 *= fInvLength

    # compute q2
    fDot1 = m01*m02 + m11*m12 + m21*m22

    fDot0 = m00*m02 + m10*m12 + m20*m22

    m02 -= fDot0*m00 + fDot1*m01
    m12 -= fDot0*m10 + fDot1*m11
    m22 -= fDot0*m20 + fDot1*m21

    fInvLength = 1.0/C.sqrt(m02*m02 + m12*m12 + m22*m22)

    m02 *= fInvLength
    m12 *= fInvLength
    m22 *= fInvLength

    return C.vertcat([C.horzcat([m00,m01,m02]),
                      C.horzcat([m10,m11,m12]),
                      C.horzcat([m20,m21,m22])])
