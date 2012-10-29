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
