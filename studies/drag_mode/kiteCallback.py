from autogen.tokiteProto import toProto
from autogen.kite_pb2 import kiteTrajectory

def callback(traj,myiter,ocp,conf,showAllPoints=False):
    kiteProtos = []
    for k in range(0,ocp.nk):
        for nicpIdx in range(0,ocp.nicp):
            if showAllPoints:
                degIdxRange = range(ocp.deg+1)
            else:
                degIdxRange = [0]
            for degIdx in degIdxRange:
                lookup = lambda name: traj.lookup(name,timestep=k,nicpIdx=nicpIdx,degIdx=degIdx)
                kiteProtos.append( toProto(lookup) )
    dkt = kiteTrajectory()
    dkt.traj.extend(list(kiteProtos))
    dkt.iteration = myiter

    dkt.messages.append("iter: "+str(myiter))
    return dkt.SerializeToString()
