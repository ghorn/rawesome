import casadi as C

def getCosLineAngle(ocp,k):
    r31 = ocp.lookup('e31',timestep=k)
    r32 = ocp.lookup('e32',timestep=k)
    r33 = ocp.lookup('e33',timestep=k)

    x = ocp.lookup('x',timestep=k)
    y = ocp.lookup('y',timestep=k)
    z = ocp.lookup('z',timestep=k)
    
#    r = ocp.lookup('r',timestep=k)
    r = C.sqrt(x*x + y*y + z*z)
    
    return (r31*x + r32*y + r33*z)/r

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
    (yaw0,pitch0,roll0) = kiteutils.getEuler(ocp, 0)
    (yawF,pitchF,rollF) = kiteutils.getEuler(ocp, -1)
    ocp.constrain(yaw0,'==',yawF)
    ocp.constrain(pitch0,'==',pitchF)
    ocp.constrain(roll0,'==',rollF)

def getOrthonormalizedDcm(ocp,k):
    m = {}
    m['e11'] = ocp.lookup('e11',timestep=k)
    m['e12'] = ocp.lookup('e12',timestep=k)
    m['e13'] = ocp.lookup('e13',timestep=k)

    m['e21'] = ocp.lookup('e21',timestep=k)
    m['e22'] = ocp.lookup('e22',timestep=k)
    m['e23'] = ocp.lookup('e23',timestep=k)

    m['e31'] = ocp.lookup('e31',timestep=k)
    m['e32'] = ocp.lookup('e32',timestep=k)
    m['e33'] = ocp.lookup('e33',timestep=k)

    return orthonormalizeDcm(m)

# dcm periodic constraints
def periodicDcm(ocp):
    dcm0 = getOrthonormalizedDcm(ocp, 0)
    dcmf = getOrthonormalizedDcm(ocp, -1)
    ocp.constrain(dcm0['e11'], '==', dcmf['e11'])
    ocp.constrain(dcm0['e22'], '==', dcmf['e22'])
    ocp.constrain(dcm0['e33'], '==', dcmf['e33'])

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

    m00 = m['e11']
    m01 = m['e12']
    m02 = m['e13']

    m10 = m['e21']
    m11 = m['e22']
    m12 = m['e23']

    m20 = m['e31']
    m21 = m['e32']
    m22 = m['e33']
    
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

    return {'e11':m00,
            'e12':m01,
            'e13':m02,
            
            'e21':m10,
            'e22':m11,
            'e23':m12,
            
            'e31':m20,
            'e32':m21,
            'e33':m22
            }
