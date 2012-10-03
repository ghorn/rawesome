from IPython.core.debugger import Tracer; debug_here = Tracer()
import casadi as C

def pendulum_model(endTimeSteps=None):
    zNames =[ "ddx"
            , "ddz"
            , "tau"
            ]
    xNames = [ "x"
             , "z"
             , "dx"
             , "dz"
             ]
    
    uNames = [ "torque"
             ]

    pNames = [ "m"
             ]
    r = 0.3
    
    uSyms = [C.ssym(name) for name in uNames]
    uVec = C.veccat( uSyms )
    u = dict(zip(uNames,uSyms))

    pSyms = [C.ssym(name) for name in pNames]
    pVec = C.veccat( pSyms )
    p = dict(zip(pNames,pSyms))

    xSyms = [C.ssym(name) for name in xNames]
    xVec = C.veccat( xSyms )
    x = dict(zip(xNames,xSyms))

    zSyms = [C.ssym(name) for name in zNames]
    zVec = C.veccat( zSyms )
    z = dict(zip(zNames, zSyms))

    stateDotDummy = C.veccat( [C.ssym(name+"DotDummy") for name in xNames] )

    scaledStateDotDummy = stateDotDummy
    
    if endTimeSteps is not None:
        endTime,nSteps = endTimeSteps
        pNames.append("endTime")
        p["endTime"] = endTime
        pVec = C.veccat([pVec,endTime])
        scaledStateDotDummy = stateDotDummy/(endTime/(nSteps-1))


    xDotDummy  = scaledStateDotDummy[0]
    zDotDummy  = scaledStateDotDummy[1]
    dxDotDummy = scaledStateDotDummy[2]
    dzDotDummy = scaledStateDotDummy[3]

    ode = [ x['dx'] - xDotDummy
          , x['dz'] - zDotDummy
          , z['ddx'] - dxDotDummy
          , z['ddz'] - dzDotDummy
          ]

    fx =  u['torque']*x['z']
    fz = -u['torque']*x['x'] + p['m']*9.8
    alg = [ p['m']*z['ddx'] + x['x']*z['tau'] - fx,
            p['m']*z['ddz'] + x['z']*z['tau'] - fz,
            x['x']*z['ddx'] + x['z']*z['ddz'] + (x['dx']*x['dx'] + x['dz']*x['dz']) ]

    c = [ x['x']*x['x'] + x['z']*x['z'] - r*r
        , x['dx']*x['x']* + x['dz']*x['z']
        ]
    
    dae = C.SXFunction( C.daeIn( x=xVec, z=zVec, p=C.veccat([uVec,pVec]), xdot=stateDotDummy ),
                        C.daeOut( alg=C.veccat(alg), ode=C.veccat(ode)))
    return (dae, {'xVec':xVec,'xDict':x,'xNames':xNames,
                  'uVec':uVec,'uNames':uNames,
                  'pVec':pVec,'pNames':pNames,
                  'zVec':zVec,'c':c})

if __name__=='__main__':
    (f,others) = rocket_model()
