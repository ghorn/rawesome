import kite_pb2
import numpy

def toKiteProto(lookup,kiteAlpha=1.0,lineAlpha=1.0):
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

    try:
        cs.rArm = lookup['rArm']
    except Exception:
        cs.rArm = 0.0
        pass
    try:
        cs.zt = lookup['zt']
    except Exception:
        cs.zt = 0
        pass

    cs.w0 = lookup('w0')
        
    cs.kiteTransparency = kiteAlpha
    cs.lineTransparency = lineAlpha

    for (attrName,lookupName) in [('CL','cL'),
                                  ('CD','cD'),
                                  ('L_over_D','L/D'),
                                  ('alpha_deg','alpha(deg)'),
                                  ('beta_deg','beta(deg)'),
                                  ('airspeed','airspeed'),
                                  ('tension','tether tension'),
                                  ('power','winch power'),
                                  ('energy','quadrature energy'),
                                  ('line_angle_deg', 'line angle (deg)'),
                                  ('r',  'r'),
                                  ('dr', 'dr'),
                                  ('ddr','ddr'),
                                  ('c','c'),
                                  ('cdot','cdot'),
                                  ('elevator_deg','elevator(deg)'),
                                  ('aileron_deg','aileron(deg)')]:
        try:
            cs.outputs.__setattr__(attrName,lookup(lookupName))
        except Exception:
            pass
    return cs
