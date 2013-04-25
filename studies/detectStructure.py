import rawe
import casadi as C
import numpy

import Numberjack as nj

if __name__=='__main__':
    print "creating model..."
    from highwind_carousel_conf import conf
    dae = rawe.models.carousel(conf)
    #dae = rawe.models.pendulum2()

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
        
    def decompose(A):
      A = C.IMatrix(A)
      C.makeSparse(A)
      rowperm = C.IVector()
      colperm = C.IVector()
      rowblock = C.IVector()
      colblock = C.IVector()
      coarse_rowblock = C.IVector()
      coarse_colblock = C.IVector()

      A.sparsity().dulmageMendelsohn 	( rowperm, colperm, rowblock,	colblock,	coarse_rowblock, coarse_colblock)

      print "rowperm: ", rowperm
      print "colperm: ", colperm
      print "recovered structure:"
      A[rowperm,colperm].printMatrix()
      print "rowblock: ", rowblock
      print "colblock: ", colblock
      print "coarse_rowblock: ", coarse_rowblock
      print "coarse_colblock: ", coarse_colblock
      return rowperm,colperm
        
        
    MA = qualifyM(jac[:,:nx])
    MZ = qualifyM(jac[:,nx:nx+nz])
    MU = qualifyM(jac[:,nx+nz:nx+nz+nup])
    MC = qualifyM(jac[:,nx+nz+nup:])
    
    
    X = nj.VarArray(nx,3)  
    Y = nj.VarArray(nx+nz,3)
    
    model = nj.Model(nj.Minimise( nj.Sum([x==1 for x in X]) ))
   # model = nj.Model(nj.Minimise( nj.Sum(X) ))
    # Contraints on A
    # 1 0 0
    # 2 2 0
    # 2 2 1
    for i in range(MA.size1()):
      for j in range(MA.size2()):
        e = MA[i,j].at(0)
        if e==0:
          pass
        elif e==1:
          model+= [ 2-Y[i]+X[j]<=2 ]
        elif e==2:
          model+=[ Y[i]>=1, X[j]<=1]
        else:
          raise Exception("?")

    # Contraints on C
    # 1 0 0
    # 2 2 0
    # 2 2 1
    for i in range(MC.size1()):
      for j in range(MC.size2()):
        e = MC[i,j].at(0)
        if e==0:
          pass
        elif e==1:
          model+= [ 2-Y[i]+X[j]<=2 ]
        elif e==2:
          model+=[ Y[i]>=1, X[j]<=1]
        else:
          raise Exception("?")

    # Contraints on Z
    # 0
    # 2
    # 2
    for i in range(MZ.size1()):
      for j in range(MZ.size2()):
        e = MZ[i,j].at(0)
        if e==0:
          pass
        elif e==1 or e==2:
          model+= [ Y[i]>=1 ]
        else:
          raise Exception("?")

    # Contraints on U
    # 1
    # 2
    # 2
    for i in range(MU.size1()):
      for j in range(MU.size2()):
        e = MU[i,j].at(0)
        if e==0 or e==1:
          pass
        elif e==2:
          model+= [ Y[i]>=1 ]
        else:
          raise Exception("?")
          
          
    # C1 and C3 should be square 
    model+= [ nj.Sum([x==0 for x in X]) == nj.Sum([y==0 for y in Y])]
    model+= [ nj.Sum([x==2 for x in X]) == nj.Sum([y==2 for y in Y])]
    
    #import SCIP 
    #import MiniSat
    import Mistral
    #import Walksat
          
    solver = Mistral.Solver(model, [X,Y])
    solver.setVerbosity(2)
    solver.solve()
    
    Xi =  [[i for (i,e) in enumerate(X) if e.get_value()==j] for j in range(3)]
    Yi =  [[i for (i,e) in enumerate(Y) if e.get_value()==j] for j in range(3)]
 
    def blockprint(A,brow,bcol):
      assert sum(brow)==A.size1()
      assert sum(bcol)==A.size2()
      bcol = numpy.cumsum(bcol)
      brow = numpy.cumsum(brow)
      s=""
      for i in range(A.size1()):
        s+= sum([i==k for k in brow])*("-" * ((A.size2()+len(bcol)-1)*2-1) + "\n")
        for j in range(A.size2()):
          s+= sum([j==k for k in bcol])*"| "
          m = A[i,j].at(0)
          s+="%d " % m
        s+= max(sum([A.size2()==k for k in bcol])-1,0)*"| " 
        s+="\n"
      s+=  max(sum([A.size1()==k for k in brow])-1,0)*("-" * ((A.size2()+len(bcol)-1)*2-1) + "\n")
      print s
      
    blockprint(C.blockcat([[i for i in range(10)] for j in range(12)]),[4,5,3],[1,2,7])
        
    

    print "MA", MA.dimString()
    blockprint(MA[sum(Yi,[]),sum(Xi,[])],[len(i) for i in Yi],[len(i) for i in Xi])
    print "MC", MC.dimString()
    blockprint(MC[sum(Yi,[]),sum(Xi,[])],[len(i) for i in Yi],[len(i) for i in Xi])
    print "MZ", MZ.dimString()
    blockprint(MZ[sum(Yi,[]),:],[len(i) for i in Yi],[MZ.size2()])
    print "MU", MU.dimString()
    blockprint(MU[sum(Yi,[]),:],[len(i) for i in Yi],[MU.size2()])
    
    print X
    print Y
    print Xi
    print Yi
    
    for i,l in enumerate(Xi):
      print "Cat #", i, [dae.xNames()[k] for k in l]
    
