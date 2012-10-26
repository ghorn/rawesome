import kite_pb2

def toKiteProto(x,u,p,zt,rArm,w0=None,zeroDelta=False):
    cs = kite_pb2.CarouselState()

    cs.kiteXyz.x = x.at(0)
    cs.kiteXyz.y = x.at(1)
    cs.kiteXyz.z = x.at(2)

    cs.kiteDcm.r11 = x.at(3)
    cs.kiteDcm.r12 = x.at(4)
    cs.kiteDcm.r13 = x.at(5)

    cs.kiteDcm.r21 = x.at(6)
    cs.kiteDcm.r22 = x.at(7)
    cs.kiteDcm.r23 = x.at(8)

    cs.kiteDcm.r31 = x.at(9)
    cs.kiteDcm.r32 = x.at(10)
    cs.kiteDcm.r33 = x.at(11)

    if zeroDelta:
        cs.delta = 0
    else:
        cs.delta = x.at(18)

    cs.rArm = rArm
    cs.zt = zt

    if w0 is not None:
        cs.w0 = w0

    return cs
