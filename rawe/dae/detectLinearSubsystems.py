import rawe
import casadi as C

def detectLinearSubsystems(dae):
    f = dae.getResidual()
    xdot = C.veccat([dae.ddt(name) for name in dae.xNames()])
    x = dae.xVec()
    z = dae.zVec()
    p = dae.pVec()
    u = dae.uVec()
    nx = x.size()
    nz = z.size()
    nup = u.size() + p.size()

    inputs = C.veccat([x, z, u, p, xdot])

    # take jacobian and find which entries are constant, zero, or nonzer
    jac = C.jacobian(f,inputs)
    
    def qualifyM(m):
      M = C.IMatrix.zeros(m.size1(),m.size2())
      for i in range(m.size1()):
        for j in range(m.size2()):
          M[i,j] = qualify(m[i,j].toScalar())
      return M
      
    def qualify(e):
      if e.isZero():
        return 0
      f = C.SXFunction([p],[e])
      f.init()
      if len(f.getFree()) > 0:
        return 2
      else:
        return 1
        
    MA = qualifyM(jac[:,:nx])
    MZ = qualifyM(jac[:,nx:nx+nz])
    MU = qualifyM(jac[:,nx+nz:nx+nz+nup])
    MC = qualifyM(jac[:,nx+nz+nup:])
    M = qualifyM(jac)

    # which equations have nonlinearity
    fi_nonlinear = set()
    fi_linear = set()
    for i in range(f.shape[0]):
        if any(M[i,:] > 1) or any(MZ[i,:] > 0):
            fi_nonlinear.add(i)
        else:
            fi_linear.add(i)

    # which variables are in the linear set
    xj_linear = set()
    for i in fi_linear:
        for j in range(nx):
            if MA[i,j] == 2 or MC[i,j] == 2:
                raise Exception('the "impossible" happend')
#            if MC[i,j] == 1 or MA[i,j] == 1:
            if MC[i,j] == 1: # poor man's C1 square constraint
                xj_linear.add(j)

    # which variables are in the nonlinear set
    xj_nonlinear = set()
    for i in fi_nonlinear:
        for j in range(nx):
            if MA[i,j] >= 1 or MC[i,j] >= 1:
                xj_nonlinear.add(j)

    x23_candidates = xj_nonlinear - xj_linear
    x1_candidates = set(range(nx))-x23_candidates

    # kick out linear variables which depend on nonlinear ones
    changed = True
    niters = 0
    while changed == True:
        niters += 1
        changed = False
        badRows = set()
        for i in fi_linear:
            if any(MZ[i,:] > 0):
                raise Exception('the "impossible" happened')
            for j in x23_candidates:
                if MA[i,j] > 0 or MC[i,j] > 0:
                    badRows.add(i)
            
            # oh shit, a tainted row, blow away everything here
            if i in badRows:
                removeUs = set()
                for j in x1_candidates:
                    if MA[i,j] > 0 or MC[i,j] > 0:
                        changed = True
                        removeUs.add(j)
                for j in removeUs:
                    x1_candidates.remove(j)
                    x23_candidates.add(j)
        for i in badRows:
            fi_nonlinear.add(i)
            fi_linear.remove(i)
    f1 = fi_linear
    f23 = fi_nonlinear
    x1 = x1_candidates
    x23 = x23_candidates
    print "finished in",niters,"iterations"
    print "LOL HERE IS WHERE WE CHECK IF C1 IS SQUARE (yolo)"

    # separate x23 into x2 and x3
    x3_candidates = set(x23)
    for i in f23:
        not_x3 = set()
        for j in x3_candidates:
            if MC[i,j] == 2 or MA[i,j] == 2:
                not_x3.add(j)
        for j in not_x3:
            x3_candidates.remove(j)
    x3 = x3_candidates
    x2 = x23 - x3

    # separate f23 into f2 and f3
    not_f3 = set()
    for i in f23:
        for j in x3:
            if MA[i,j] > 0 or MC[i,j] > 0:
                not_f3.add(i)
    f2 = f23 - not_f3
    f3 = f23 - f2

    print "LOL HERE IS WHERE WE CHECK IF C3 IS SQUARE (yolo)"
    print "x1",[dae.xNames()[j] for j in x1]
    print "x2",[dae.xNames()[j] for j in x2]
    print "x3",[dae.xNames()[j] for j in x3]
    print "f1",f1
    print "f2",f2
    print "f3",f3
