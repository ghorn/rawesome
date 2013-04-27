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

    try:
        cs.w0 = lookup('w0')
    except Exception:
        pass
        
    cs.kiteTransparency = kiteAlpha
    cs.lineTransparency = lineAlpha

    # set forces torques
    for (attrName,lookupName) in [('fx','aero_fx'),
                                  ('fy','aero_fy'),
                                  ('fz','aero_fz'),
                                  ('mx','aero_mx'),
                                  ('my','aero_my'),
                                  ('mz','aero_mz')]:
        cs.aero.forcesTorques.__setattr__(attrName,lookup(lookupName))
    for (attrName,lookupName) in [('alpha_deg','alpha_deg'),
                                  ('beta_deg','beta_deg'),
                                  ('airspeed','airspeed'),
                                  ('CL','cL'),
                                  ('CD','cD'),
                                  ('L_over_D','L_over_D'),
                                  ('fLift','fLift'),
                                  ('fDrag','fDrag')]:
        try:
            cs.aero.__setattr__(attrName,lookup(lookupName))
        except Exception:
            pass
    for (attrName,lookupName) in [('tension','tether_tension'),
                                  ('winch_power','winch_power'),
                                  ('prop_power','prop_power'),
                                  ('energy','quadrature_energy'),
                                  ('line_angle_deg', 'line_angle_deg'),
                                  ('r',  'r'),
                                  ('dr', 'dr'),
                                  ('ddr','ddr'),
                                  ('c','c'),
                                  ('cdot','cdot'),
                                  ('elevator_deg','elevator_deg'),
                                  ('prop_drag','prop_drag'),
                                  ('aileron_deg','aileron_deg')]:
        try:
            cs.outputs.__setattr__(attrName,lookup(lookupName))
        except Exception:
            pass
    return cs
