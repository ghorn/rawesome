import rawe
import casadi as C

if __name__=='__main__':
    print "creating model..."
    from highwind_carousel_conf import conf
    dae = rawe.models.carousel(conf)

    print 'x',dae.xNames()
    print 'z',dae.zNames()
    print 'u',dae.uNames()
    print 'p',dae.pNames()

    f = dae.getResidual()
    xdot = C.veccat([dae.ddt(name) for name in dae.xNames()])
    x = dae.xVec()
    z = dae.zVec()
    p = dae.pVec()
    u = dae.uVec()
    nx = x.size()
    nz = z.size()
    nup = u.size() + p.size()

    print f.shape
    print 
    inputs = C.veccat([x, z, u, p, xdot])
    xIdx = range(0,nx)
    zIdx = range(nx,nx+nz)
    upIdx = range(nx+nz,nx+nz+nup)
    xDotIdx = range(nx+nz+nup,nx+nz+nup+nx)

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
            if MC[i,j] == 1:
                xj_linear.add(j)

    # which variables are in the nonlinear set
    xj_nonlinear = set()
    for i in fi_nonlinear:
        for j in range(nx):
            if MA[i,j] >= 1 or MC[i,j] >= 1:
                xj_nonlinear.add(j)

    print fi_nonlinear
    print fi_linear
    print [dae.xNames()[j] for j in xj_nonlinear]
    print [dae.xNames()[j] for j in xj_linear]
    print "="*80
    x23_candidates = xj_nonlinear - xj_linear
    x1_candidates = set(range(nx))-x23_candidates
    print "x23_c",[dae.xNames()[j] for j in x23_candidates]
    print "x1_c",[dae.xNames()[j] for j in x1_candidates]

    for k in fi_linear:
        print f[k]

    # kick out linear variables which depend on nonlinear ones
    changed = True
    niters = 0
    while changed == True:
        niters += 1
        changed = False
        for i in fi_linear:
            badColumn = False
            if any(MZ[i,:] > 0):
                raise Exception('the "impossible" happened')
            for j in x23_candidates:
                if MA[i,j] > 0 or MC[i,j] > 0:
                    badColumn = True
            
            # oh shit, a tainted row, blow away everything here
            if badColumn:
                print 'got a bad column, yo'
                removeUs = set()
                for j in x1_candidates:
                    if MA[i,j] > 0 or MC[i,j] > 0:
                        changed = True
                        removeUs.add(j)
                for j in removeUs:
                    x1_candidates.remove(j)
                    x23_candidates.add(j)
    print "finished in",niters,"iterations"
    print "-"*80
    print "x23_c",[dae.xNames()[j] for j in x23_candidates]
    print "x1_c",[dae.xNames()[j] for j in x1_candidates]
