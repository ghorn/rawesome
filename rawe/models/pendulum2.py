from dae import Dae
import casadi as C

def pendulumModel(nSteps=None):
    dae = Dae()
    dae.addZ( [ "tau"
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
    
    dae.stateDotDummy = C.veccat( [C.ssym(name+"DotDummy") for name in dae._xNames] )

    scaledStateDotDummy = dae.stateDotDummy
    
    if nSteps is not None:
        endTime = dae.addP('endTime')
        scaledStateDotDummy = dae.stateDotDummy/(endTime/(nSteps-1))

    xDotDummy  = scaledStateDotDummy[0]
    zDotDummy  = scaledStateDotDummy[1]
    dxDotDummy = scaledStateDotDummy[2]
    dzDotDummy = scaledStateDotDummy[3]
    
    ddx = dxDotDummy
    ddz = dzDotDummy

    ode = [ dae['dx'] - xDotDummy
          , dae['dz'] - zDotDummy
#          , dae['ddx'] - dxDotDummy
#          , dae['ddz'] - dzDotDummy
          ]

    fx =  dae['torque']*dae['z']
    fz = -dae['torque']*dae['x'] + dae['m']*9.8
    
    alg = [ dae['m']*ddx + dae['x']*dae['tau'] - fx
          , dae['m']*ddz + dae['z']*dae['tau'] - fz
          , dae['x']*ddx + dae['z']*ddz + (dae['dx']*dae['dx'] + dae['dz']*dae['dz']) ]

    dae['c']    = dae['x']*dae['x'] + dae['z']*dae['z'] - r*r
    dae['cdot'] = dae['dx']*dae['x'] + dae['dz']*dae['z']

    dae.setAlgRes( alg )
    dae.setOdeRes( ode )

    return dae

if __name__=='__main__':
    pendulumModel()
    pendulumModel(nSteps=10)
