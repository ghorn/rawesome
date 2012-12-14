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
    
    ddx = dae.ddt('dx')
    ddz = dae.ddt('dz')

    ode = [ dae['dx'] - dae.ddt('x')
          , dae['dz'] - dae.ddt('z')
#          , dae['ddx'] - dae.ddt('dx')
#          , dae['ddz'] - dae.ddt('dz')
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
