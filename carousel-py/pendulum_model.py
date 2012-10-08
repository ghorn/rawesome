from dae import Dae
import casadi as C

def pendulum_model(nSteps=None):
    dae = Dae()
    dae.addZ( [ "ddx"
              , "ddz"
              , "tau"
              ] )
    dae.addX( [ "x"
              , "z"
              , "dx"
              , "dz"
              ] )
    dae.addU( [ "torque"
              ] )
    dae.addP( [ "m"
              ] )

    r = 0.3
    
    stateDotDummy = C.veccat( [C.ssym(name+"DotDummy") for name in dae._xNames] )

    scaledStateDotDummy = stateDotDummy
    
    if nSteps is not None:
        endTime = dae.addP('endTime')
        scaledStateDotDummy = stateDotDummy/(endTime/(nSteps-1))

    xDotDummy  = scaledStateDotDummy[0]
    zDotDummy  = scaledStateDotDummy[1]
    dxDotDummy = scaledStateDotDummy[2]
    dzDotDummy = scaledStateDotDummy[3]

    ode = [ dae.x('dx') - xDotDummy
          , dae.x('dz') - zDotDummy
          , dae.z('ddx') - dxDotDummy
          , dae.z('ddz') - dzDotDummy
          ]

    fx =  dae.u('torque')*dae.x('z')
    fz = -dae.u('torque')*dae.x('x') + dae.p('m')*9.8
    alg = [ dae.p('m')*dae.z('ddx') + dae.x('x')*dae.z('tau') - fx
          , dae.p('m')*dae.z('ddz') + dae.x('z')*dae.z('tau') - fz
          , dae.x('x')*dae.z('ddx') + dae.x('z')*dae.z('ddz') + (dae.x('dx')*dae.x('dx') + dae.x('dz')*dae.x('dz')) ]

    c = [ dae.x('x')*dae.x('x') + dae.x('z')*dae.x('z') - r*r
        , dae.x('dx')*dae.x('x')* + dae.x('dz')*dae.x('z')
        ]
    dae.addOutput('invariants',C.veccat(c))

    dae.sxfun = C.SXFunction( C.daeIn( x=dae.xVec(), z=dae.zVec(), p=C.veccat([dae.uVec(),dae.pVec()]), xdot=stateDotDummy ),
                              C.daeOut( alg=C.veccat(alg), ode=C.veccat(ode)) )

    return dae

if __name__=='__main__':
    pendulum_model()
    pendulum_model(nSteps=10)
