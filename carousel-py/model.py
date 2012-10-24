import casadi as C
from dae import Dae

def forcesTorques(dae,conf):
    rho = conf['env']['rho']
    rA = conf['carousel']['rArm']
    alpha0 = conf['aero']['alpha0deg']*C.pi/180 
    
    #ROLL DAMPING
    rD = conf['aero']['rD']
    pD = conf['aero']['pD']
    yD = conf['aero']['yD']
    
    #WIND-TUNNEL PARAMETERS
    #Lift (report p. 67)
    cLA = conf['aero']['cLA']
    
    cLe = conf['aero']['cLe']

    cL0 = conf['aero']['cL0']
    
    #Drag (report p. 70)
    cDA = conf['aero']['cDA']
    cDA2 = conf['aero']['cDA2']
    cDB2 = conf['aero']['cDB2']

    cD0 = conf['aero']['cD0']
    
    #Roll (report p. 72)
    cRB  = conf['aero']['cRB']
    cRAB = conf['aero']['cRAB']
    cRr  = conf['aero']['cRr']
    
    #Pitch (report p. 74)
    cPA = conf['aero']['cPA']
    cPe = conf['aero']['cPe']
    
    cP0 = conf['aero']['cP0']
    
    #Yaw (report p. 76)
    cYB = conf['aero']['cYB']
    cYAB = conf['aero']['cYAB']

    #TAIL LENGTH
    lT = conf['kite']['lT']
    
    AR = conf['kite']['AR']
    area = conf['kite']['area']
    
    span = C.sqrt(area*AR)
    chord = C.sqrt(area/AR)
    
    ###########     model integ ###################
    x =   dae.x('x')
    y =   dae.x('y')
    z =   dae.x('z')

    e11 = dae.x('e11')
    e12 = dae.x('e12')
    e13 = dae.x('e13')

    e21 = dae.x('e21')
    e22 = dae.x('e22')
    e23 = dae.x('e23')

    e31 = dae.x('e31')
    e32 = dae.x('e32')
    e33 = dae.x('e33')
                   
    dx  =  dae.x('dx')
    dy  =  dae.x('dy')
    dz  =  dae.x('dz')

    w1  =  dae.x('w1')
    w2  =  dae.x('w2')
    w3  =  dae.x('w3')

    delta = dae.x('delta')
    ddelta = dae.x('ddelta')

    u1 = dae.u('aileron')
    u2 = dae.u('elevator')
    
    # wind
    z0 = conf['wind shear']['z0']
    zt_roughness = conf['wind shear']['zt_roughness']
    zsat = 0.5*(z+C.sqrt(z*z))
    wind_x = dae.p('w0')*C.log((zsat+zt_roughness+2)/zt_roughness)/C.log(z0/zt_roughness)
    dae.addOutput('wind at altitude', wind_x)
    dae.addOutput('w0', dae.p('w0'))

    dp_carousel_frame = C.veccat( [ dx - ddelta*y
                                  , dy + ddelta*(rA + x)
                                  , dz
                                  ]) - C.veccat([C.cos(delta)*wind_x,C.sin(delta)*wind_x,0])
    R_c2b = C.veccat( [e11, e12, e13,
                       e21, e22, e23,
                       e31, e32, e33] ).reshape((3,3))
    dpE = C.mul( R_c2b, dp_carousel_frame )
    
    ##### more model_integ ###########
    # EFFECTIVE WIND IN THE KITE`S SYSTEM :
    # ###############################
    
    #Airfoil speed in carousel frame
    we1 = dp_carousel_frame[0]
    we2 = dp_carousel_frame[1]
    we3 = dp_carousel_frame[2]
    
    vKite2 = C.mul(dp_carousel_frame.trans(), dp_carousel_frame) #Airfoil speed^2 
    vKite = C.sqrt(vKite2) #Airfoil speed
    dae.addOutput('airspeed', vKite)
    
    # CALCULATION OF THE FORCES :
    # ###############################
    #
    #   FORCE ARE COMPUTED IN THE CAROUSEL FRAME @>@>@>
    
    #Aero coeff.
    
            
    # LIFT DIRECTION VECTOR
    # ############
    
    # Relative wind speed in Airfoil's referential 'E'
    wE1 = dpE[0]
    wE2 = dpE[1]
    wE3 = dpE[2]
    
    # Airfoil's transversal axis in carousel referential 'e'
    eTe1 = e21
    eTe2 = e22
    eTe3 = e23
    
    
    # Lift axis ** Normed to we @>@> **
    eLe1 = - eTe2*we3 + eTe3*we2
    eLe2 = - eTe3*we1 + eTe1*we3
    eLe3 = - eTe1*we2 + eTe2*we1
       
    
    # AERODYNAMIC COEEFICIENTS
    # #################
    vT1 =          wE1
    vT2 = -lT*w3 + wE2
    vT3 =  lT*w2 + wE3
    
    
    alpha = alpha0-wE3/wE1
    
    #NOTE: beta & alphaTail are compensated for the tail motion induced by
    #omega @>@>
    betaTail = vT2/C.sqrt(vT1*vT1 + vT3*vT3)
    beta = wE2/C.sqrt(wE1*wE1 + wE3*wE3)
    alphaTail = alpha0-vT3/vT1
    
    dae.addOutput('alpha(deg)', alpha*180/C.pi)
    dae.addOutput('alphaTail(deg)', alphaTail*180/C.pi)
    dae.addOutput('beta(deg)', beta*180/C.pi)
    dae.addOutput('betaTail(deg)', betaTail*180/C.pi)

    # cL = cLA*alpha + cLe*u2   + cL0
    # cD = cDA*alpha + cDA2*alpha*alpha + cDB2*beta*beta + cDe*u2 + cDr*u1 + cD0
    # cR = -rD*w1 + cRB*beta + cRAB*alphaTail*beta + cRr*u1
    # cP = cPA*alphaTail + cPe*u2  + cP0
    # cY = cYB*beta + cYAB*alphaTail*beta
    
    cL = cLA*alpha + cLe*u2 + cL0
    cD = cDA*alpha + cDA2*alpha*alpha + cDB2*beta*beta + cD0
    cR = -rD*w1 + cRB*betaTail + cRr*u1 + cRAB*alphaTail*betaTail
    cP = -pD*w2 + cPA*alphaTail + cPe*u2 + cP0
    cY = -yD*w3 + cYB*betaTail + cYAB*alphaTail*betaTail

    dae.addOutput('cL', cL)
    dae.addOutput('cD', cD)
    dae.addOutput('L/D', cL/cD)
    
    # LIFT :
    # ###############################
    fL1 =  rho*cL*eLe1*vKite/2.0
    fL2 =  rho*cL*eLe2*vKite/2.0
    fL3 =  rho*cL*eLe3*vKite/2.0
    
    # DRAG :
    # #############################
    fD1 = -rho*vKite*cD*we1/2.0
    fD2 = -rho*vKite*cD*we2/2.0 
    fD3 = -rho*vKite*cD*we3/2.0 
    
    
    
    # FORCES (AERO)
    # ###############################
    f1 = fL1 + fD1
    f2 = fL2 + fD2
    f3 = fL3 + fD3
    
    #f = f-fT
       
    # TORQUES (AERO)
    # ###############################
     
    t1 =  0.5*rho*vKite2*span*cR
    t2 =  0.5*rho*vKite2*chord*cP
    t3 =  0.5*rho*vKite2*span*cY
    return (f1, f2, f3, t1, t2, t3)
    

def modelInteg(dae, conf):
    #  PARAMETERS OF THE KITE :
    #  ##############
    m =  conf['kite']['mass'] #  mass of the kite               #  [ kg    ]
                 
    #   PHYSICAL CONSTANTS :
    #  ##############
    g = conf['env']['g'] #  gravitational constant         #  [ m /s^2]
    
    #  PARAMETERS OF THE CABLE :
    #  ##############
     
    #INERTIA MATRIX (Kurt's direct measurements)
    j1 =  conf['kite']['j1']
    j31 = conf['kite']['j31']
    j2 =  conf['kite']['j2']
    j3 =  conf['kite']['j3']
    
    #Carousel Friction & inertia
    jCarousel = conf['carousel']['jCarousel']
    cfric = conf['carousel']['cfric']

    zt = conf['kite']['zt']
    rA = conf['carousel']['rArm']

    ###########     model integ ###################
    x =   dae.x('x')
    y =   dae.x('y')
    z =   dae.x('z')

    e11 = dae.x('e11')
    e12 = dae.x('e12')
    e13 = dae.x('e13')

    e21 = dae.x('e21')
    e22 = dae.x('e22')
    e23 = dae.x('e23')

    e31 = dae.x('e31')
    e32 = dae.x('e32')
    e33 = dae.x('e33')
                   
    dx  =  dae.x('dx')
    dy  =  dae.x('dy')
    dz  =  dae.x('dz')

    w1  =  dae.x('w1')
    w2  =  dae.x('w2')
    w3  =  dae.x('w3')

    ddelta = dae.x('ddelta')

    r = dae.x('r')
    dr = dae.x('dr')
    
    ddr = dae.u('ddr')
    
    tc = dae.u('tc') #Carousel motor torque

    (f1, f2, f3, t1, t2, t3) = forcesTorques(dae,conf)

    # ATTITUDE DYNAMICS
    # #############################
    
    # DynFile          # Call DAE
#    mm :: Matrix (Expr Double)
    mm = C.SXMatrix(8, 8)
    mm[0,0] = jCarousel + m*rA*rA + m*x*x + m*y*y + 2*m*rA*x 
    mm[0,1] = -m*y 
    mm[0,2] = m*(rA + x) 
    mm[0,3] = 0 
    mm[0,4] = 0 
    mm[0,5] = 0 
    mm[0,6] = 0 
    mm[0,7] = 0
    
    mm[1,0] = -m*y 
    mm[1,1] = m 
    mm[1,2] = 0 
    mm[1,3] = 0 
    mm[1,4] = 0 
    mm[1,5] = 0 
    mm[1,6] = 0 
    mm[1,7] = x + zt*e31
    
    mm[2,0] = m*(rA + x) 
    mm[2,1] = 0 
    mm[2,2] = m 
    mm[2,3] = 0 
    mm[2,4] = 0 
    mm[2,5] = 0 
    mm[2,6] = 0 
    mm[2,7] = y + zt*e32
    
    mm[3,0] = 0 
    mm[3,1] = 0 
    mm[3,2] = 0 
    mm[3,3] = m 
    mm[3,4] = 0 
    mm[3,5] = 0 
    mm[3,6] = 0 
    mm[3,7] = z + zt*e33
    
    mm[4,0] = 0 
    mm[4,1] = 0 
    mm[4,2] = 0 
    mm[4,3] = 0 
    mm[4,4] = j1 
    mm[4,5] = 0 
    mm[4,6] = j31 
    mm[4,7] = -zt*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33)
    
    mm[5,0] = 0 
    mm[5,1] = 0 
    mm[5,2] = 0 
    mm[5,3] = 0 
    mm[5,4] = 0 
    mm[5,5] = j2 
    mm[5,6] = 0 
    mm[5,7] = zt*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33)
    
    mm[6,0] = 0 
    mm[6,1] = 0 
    mm[6,2] = 0 
    mm[6,3] = 0 
    mm[6,4] = j31 
    mm[6,5] = 0 
    mm[6,6] = j3 
    mm[6,7] = 0
    
    mm[7,0] = -zt*(e11*e23*x - e13*e21*x + e12*e23*y - e13*e22*y + zt*e11*e23*e31 - zt*e13*e21*e31 + zt*e12*e23*e32 - zt*e13*e22*e32) 
    mm[7,1] = x + zt*e31 
    mm[7,2] = y + zt*e32 
    mm[7,3] = z + zt*e33 
    mm[7,4] = -zt*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33) 
    mm[7,5] = zt*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) 
    mm[7,6] = 0 
    mm[7,7] = 0
    
    zt2 = zt*zt
    rhs = C.veccat(
          [ tc - cfric*ddelta - f1*y + f2*(rA + x) + dy*m*(dx - 2*ddelta*y) - dx*m*(dy + 2*ddelta*rA + 2*ddelta*x) 
          , f1 + ddelta*m*(dy + ddelta*rA + ddelta*x) + ddelta*dy*m 
          , f2 - ddelta*m*(dx - ddelta*y) - ddelta*dx*m 
          , f3 - g*m 
          , t1 - w2*(j3*w3 + j31*w1) + j2*w2*w3 
          , t2 + w1*(j3*w3 + j31*w1) - w3*(j1*w1 + j31*w3) 
          , t3 + w2*(j1*w1 + j31*w3) - j2*w1*w2
          , ddr*r-(zt*w1*(e11*x+e12*y+e13*z+zt*e11*e31+zt*e12*e32+zt*e13*e33)+zt*w2*(e21*x+e22*y+e23*z+zt*e21*e31+zt*e22*e32+zt*e23*e33))*(w3-ddelta*e33)-dx*(dx-zt*e21*(w1-ddelta*e13)+zt*e11*(w2-ddelta*e23))-dy*(dy-zt*e22*(w1-ddelta*e13)+zt*e12*(w2-ddelta*e23))-dz*(dz-zt*e23*(w1-ddelta*e13)+zt*e13*(w2-ddelta*e23))+dr*dr+(w1-ddelta*e13)*(e21*(zt*dx-zt2*e21*(w1-ddelta*e13)+zt2*e11*(w2-ddelta*e23))+e22*(zt*dy-zt2*e22*(w1-ddelta*e13)+zt2*e12*(w2-ddelta*e23))+zt*e23*(dz+zt*e13*w2-zt*e23*w1)+zt*e33*(w1*z+zt*e33*w1+ddelta*e11*x+ddelta*e12*y+zt*ddelta*e11*e31+zt*ddelta*e12*e32)+zt*e31*(x+zt*e31)*(w1-ddelta*e13)+zt*e32*(y+zt*e32)*(w1-ddelta*e13))-(w2-ddelta*e23)*(e11*(zt*dx-zt2*e21*(w1-ddelta*e13)+zt2*e11*(w2-ddelta*e23))+e12*(zt*dy-zt2*e22*(w1-ddelta*e13)+zt2*e12*(w2-ddelta*e23))+zt*e13*(dz+zt*e13*w2-zt*e23*w1)-zt*e33*(w2*z+zt*e33*w2+ddelta*e21*x+ddelta*e22*y+zt*ddelta*e21*e31+zt*ddelta*e22*e32)-zt*e31*(x+zt*e31)*(w2-ddelta*e23)-zt*e32*(y+zt*e32)*(w2-ddelta*e23))
          ] )
 
    dRexp = C.SXMatrix(3,3)

    dRexp[0,0] = e21*(w3 - ddelta*e33) - e31*(w2 - ddelta*e23) 
    dRexp[0,1] = e31*(w1 - ddelta*e13) - e11*(w3 - ddelta*e33) 
    dRexp[0,2] = e11*(w2 - ddelta*e23) - e21*(w1 - ddelta*e13) 

    dRexp[1,0] = e22*(w3 - ddelta*e33) - e32*(w2 - ddelta*e23) 
    dRexp[1,1] = e32*(w1 - ddelta*e13) - e12*(w3 - ddelta*e33) 
    dRexp[1,2] = e12*(w2 - ddelta*e23) - e22*(w1 - ddelta*e13) 

    dRexp[2,0] = e23*w3 - e33*w2 
    dRexp[2,1] = e33*w1 - e13*w3 
    dRexp[2,2] = e13*w2 - e23*w1
    
    c =(x + zt*e31)**2/2 + (y + zt*e32)**2/2 + (z + zt*e33)**2/2 - r**2/2
    
    cdot =dx*(x + zt*e31) + dy*(y + zt*e32) + dz*(z + zt*e33) + zt*(w2 - ddelta*e23)*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) - zt*(w1 - ddelta*e13)*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33) - r*dr
#  cddot = (zt*w1*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) + zt*w2*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33))*(w3 - ddelta*e33) + dx*(dx + zt*e11*w2 - zt*e21*w1 - zt*ddelta*e11*e23 + zt*ddelta*e13*e21) + dy*(dy + zt*e12*w2 - zt*e22*w1 - zt*ddelta*e12*e23 + zt*ddelta*e13*e22) + dz*(dz + zt*e13*w2 - zt*e23*w1) + ddx*(x + zt*e31) + ddy*(y + zt*e32) + ddz*(z + zt*e33) - (w1 - ddelta*e13)*(e21*(zt*dx - zt**2*e21*(w1 - ddelta*e13) + zt**2*e11*(w2 - ddelta*e23)) + e22*(zt*dy - zt**2*e22*(w1 - ddelta*e13) + zt**2*e12*(w2 - ddelta*e23)) + zt*e33*(z*w1 + ddelta*e11*x + ddelta*e12*y + zt*e33*w1 + zt*ddelta*e11*e31 + zt*ddelta*e12*e32) + zt*e23*(dz + zt*e13*w2 - zt*e23*w1) + zt*e31*(w1 - ddelta*e13)*(x + zt*e31) + zt*e32*(w1 - ddelta*e13)*(y + zt*e32)) + (w2 - ddelta*e23)*(e11*(zt*dx - zt**2*e21*(w1 - ddelta*e13) + zt**2*e11*(w2 - ddelta*e23)) + e12*(zt*dy - zt**2*e22*(w1 - ddelta*e13) + zt**2*e12*(w2 - ddelta*e23)) - zt*e33*(z*w2 + ddelta*e21*x + ddelta*e22*y + zt*e33*w2 + zt*ddelta*e21*e31 + zt*ddelta*e22*e32) + zt*e13*(dz + zt*e13*w2 - zt*e23*w1) - zt*e31*(w2 - ddelta*e23)*(x + zt*e31) - zt*e32*(w2 - ddelta*e23)*(y + zt*e32)) + zt*(dw2 - dddelta*e23)*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) - zt*(dw1 - dddelta*e13)*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33) - zt*dddelta*(e11*e23*x - e13*e21*x + e12*e23*y - e13*e22*y + zt*e11*e23*e31 - zt*e13*e21*e31 + zt*e12*e23*e32 - zt*e13*e22*e32)
#      where
#        dw1 = dw @> 0
#        dw2 = dw @> 1
#        {-
#        dw3 = dw @> 2
#        -}
#        ddx = ddX @> 0
#        ddy = ddX @> 1
#        ddz = ddX @> 2
#        dddelta = dddelta' @> 0

    dae.addOutput('c',  c)
    dae.addOutput('cdot', cdot)
    return (mm, rhs, dRexp, c, cdot)
        
def model(conf,nSteps=None,extraParams=[]):
    dae = Dae()
    for ep in extraParams:
        dae.addP(ep)
        
    dae.addZ( [ "dddelta"
              , "ddx"
              , "ddy"
              , "ddz"
              , "dw1"
              , "dw2"
              , "dw3"
              , "nu"
              ] )
    dae.addX( [ "x"   # state 0
              , "y"   # state 1
              , "z"   # state 2
              , "e11" # state 3
              , "e12" # state 4
              , "e13" # state 5
              , "e21" # state 6
              , "e22" # state 7
              , "e23" # state 8
              , "e31" # state 9
              , "e32" # state 10
              , "e33" # state 11
              , "dx"  # state 12
              , "dy"  # state 13
              , "dz"  # state 14
              , "w1"  # state 15
              , "w2"  # state 16
              , "w3"  # state 17
              , "delta" # state 18
              , "ddelta" # state 19
              , "r" # state 20
              , "dr" # state 21
              , "energy" # state 22
              ] )
    dae.addU( [ "tc"
              , "aileron"
              , "elevator"
              , 'ddr'
              ] )
    dae.addP( ['w0'] )
    
    dae.addOutput('r', dae.x('r'))
    dae.addOutput('dr', dae.x('dr'))
    dae.addOutput('RPM', dae.x('ddelta')*60/(2*C.pi))
    dae.addOutput('aileron(deg)', dae.u('aileron')*180/C.pi)
    dae.addOutput('elevator(deg)', dae.u('elevator')*180/C.pi)
    
    dae.addOutput('motor torque', dae.u('tc'))
    dae.addOutput('motor power', dae.u('tc')*dae.x('ddelta'))

    dae.addOutput('winch force', dae.x('r')*dae.z('nu'))
    dae.addOutput('winch power', dae.x('r')*dae.x('dr')*dae.z('nu'))
    
    (massMatrix, rhs, dRexp, c, cdot) = modelInteg(dae, conf)

    
    ode = C.veccat( [ C.veccat(dae.x(['dx','dy','dz']))
                    , dRexp.trans().reshape([9,1])
                    , C.veccat(dae.z(['ddx','ddy','ddz']))
                    , C.veccat(dae.z(['dw1','dw2','dw3']))
                    , C.veccat([dae.x('ddelta'), dae.z('dddelta')])
                    , dae.x('dr')
                    , dae.u('ddr')
                    , dae.output('winch power') + dae.output('motor power')
                    ] )

    dae.stateDotDummy = C.veccat( [C.ssym(name+"DotDummy") for name in dae._xNames] )
    scaledStateDotDummy = dae.stateDotDummy
    
    if nSteps is not None:
        endTime = dae.addP('endTime')
        scaledStateDotDummy = dae.stateDotDummy/(endTime/(nSteps-1))

    dae.setOdeRes( ode - scaledStateDotDummy )
    dae.setAlgRes( C.mul(massMatrix, dae.zVec()) - rhs )
    
    return dae

if __name__=='__main__':
    (f,others) = model()
