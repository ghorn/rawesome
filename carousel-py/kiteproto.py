import kite_pb2

def toKiteProto(x,u,p=None):
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

    cs.delta = x.at(18)
    cs.ddelta = x.at(19)

    cs.tc = u.at(0)
    cs.u1 = u.at(1)
    cs.u2 = u.at(2)
    cs.rArm = 1.085
    cs.zt = -0.01

    if u.size()>=4:
        cs.wind_x = u.at(3)
    return cs
