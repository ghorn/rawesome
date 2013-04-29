import kiteproto
import kite_pb2

def showAllPoints(traj,myiter,ocp,conf):
    return normalCallback(traj,myiter,ocp,conf,showAllPoints=True)

def normalCallback(traj,myiter,ocp,conf,showAllPoints=False):
    kiteProtos = []
    for k in range(0,ocp.nk):
        for nicpIdx in range(0,ocp.nicp):
            if showAllPoints:
                degIdxRange = range(ocp.deg+1)
            else:
                degIdxRange = [0]
            for degIdx in degIdxRange:
                lookup = lambda name: traj.lookup(name,timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)
                kiteProtos.append( kiteproto.toKiteProto(lookup,
                                                         lineAlpha=0.2) )
    mc = kite_pb2.MultiCarousel()
    mc.horizon.extend(list(kiteProtos))

    mc.messages.append("w0: "+str(traj.lookup('w0')))
    mc.messages.append("iter: "+str(myiter))
    mc.messages.append("endTime: "+str(traj.lookup('endTime')))
    mc.messages.append("average power: "+str(traj.lookup('quadrature_energy',timestep=-1)/traj.lookup('endTime'))+" W")
    
    return mc.SerializeToString()
