import rawe
import casadi as C
import numpy

if __name__=='__main__':
    print "creating model..."
    from highwind_carousel_conf import conf
    dae = rawe.models.carousel(conf)

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
    xIdx = range(0,nx)
    zIdx = range(nx,nx+nz)
    upIdx = range(nx+nz,nx+nz+nup)
    xDotIdx = range(nx+nz+nup,nx+nz+nup+nx)

    # take jacobian and find which entries are constant, zero, or nonzer
    jac = C.jacobian(f,inputs)

    nf = f.size()

    zero     = numpy.zeros((nf, inputs.size()))
    nonzero  = numpy.zeros((nf, inputs.size()))
    constant = numpy.zeros((nf, inputs.size()))
    for i in range(nf):
        for j in range(inputs.size()):
            constant[i,j] = jac[i,j].toScalar().isConstant()
            zero[i,j]     = jac[i,j].toScalar().isZero()
            nonzero[i,j]  = not zero[i,j]
#            print (i,j), "const:",const, "zero:",zero, str(jac[i,j])[:20]

    finfo = []
    for i in range(nf):
        linearList = []
        nonlinearList = []
        nonzeroList = []
        for j in range(inputs.size()):
            if nonzero[i,j]:
                nonzeroList.append(j)
                if constant[i,j]:
                    linearList.append(j)
                else:
                    nonlinearList.append(j)
        finfo.append({'linear':linearList, 'nonlinear':nonlinearList,'nonzero':nonzeroList})

    ####### find C1*x1d = A1*x1 + B1*u #####
    # get f which are linear
    f1maybe = []
    for i in range(nf):
        fi = finfo[i]
        
        l = fi['linear']
        nl = fi['nonlinear']
        nz = fi['nonzero']
        # no nonlinear fields
        if len(nl) > 0:
            continue
        # no z
        if any([zk in fi['nonzero'] for zk in zIdx]):
            continue
        
        # ok, this could be an f1
        f1maybe.append(i)
    
    numlin = len(f1maybe)
