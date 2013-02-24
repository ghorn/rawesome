import kite_pb2
import numpy

def toKiteProto(lookup,zt,rArm,w0=None,kiteAlpha=1.0,lineAlpha=1.0):
    cs = kite_pb2.CarouselState()

    cs.kiteXyz.x = lookup('x')
    cs.kiteXyz.y = lookup('y')
    cs.kiteXyz.z = lookup('z')

    cs.kiteDcm.r11 = lookup('e11')
    cs.kiteDcm.r12 = lookup('e12')
    cs.kiteDcm.r13 = lookup('e13')

    cs.kiteDcm.r21 = lookup('e21')
    cs.kiteDcm.r22 = lookup('e22')
    cs.kiteDcm.r23 = lookup('e23')

    cs.kiteDcm.r31 = lookup('e31')
    cs.kiteDcm.r32 = lookup('e32')
    cs.kiteDcm.r33 = lookup('e33')

    try:
        cs.delta = lookup('delta')
    except:
        cs.delta = 0
        pass

    cs.rArm = rArm
    cs.zt = zt

    if w0 is not None:
        cs.w0 = w0
        
    cs.kiteTransparency = kiteAlpha
    cs.lineTransparency = lineAlpha

#    cs.outputs.CL = lookup('cL')
#    cs.outputs.CD = lookup('cD')
#    cs.outputs.L_over_D = lookup('L/D')
#    cs.outputs.alpha_deg = lookup('alpha(deg)')
#    cs.outputs.beta_deg = lookup('beta(deg)')
#    cs.outputs.airspeed = lookup('airspeed')
#    # the next two aren't defined at degIdx==0
#    try:
#        cs.outputs.tension = lookup('tether tension')
#    except:
#        pass
#    try:
#        cs.outputs.power = lookup('winch power')
#    except:
#        pass
#    # this isn't defined for the homotopy
#    try:
#        cs.outputs.energy = lookup('quadrature energy')
#    except:
#        pass
#    cs.outputs.line_angle_deg = lookup('line angle (deg)')
#    cs.outputs.r   = lookup('r')
#    cs.outputs.dr  = lookup('dr')
#    # controls are only defined at the beginning of the interval
#    try:
#        cs.outputs.ddr = lookup('ddr')
#    except:
#        pass

    return cs
