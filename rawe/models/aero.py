import casadi as C

def aeroForcesTorques(dae, conf, we, wE, (w1,w2,w3), (eTe1, eTe2, eTe3), (aileron,elevator)):
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

    sref = conf['kite']['sref']
    bref = conf['kite']['bref']
    cref = conf['kite']['cref']
    
    ##### more model_integ ###########
    # EFFECTIVE WIND IN THE KITE`S SYSTEM :
    # ###############################
    
    #Airfoil speed in carousel frame
    we1 = we[0]
    we2 = we[1]
    we3 = we[2]
    
    vKite2 = C.mul(we.trans(), we) #Airfoil speed^2 
    vKite = C.sqrt(vKite2) #Airfoil speed
    dae['airspeed'] = vKite
    
    # CALCULATION OF THE FORCES :
    # ###############################
    #
    #   FORCE ARE COMPUTED IN THE CAROUSEL FRAME @>@>@>
    
            
    # LIFT DIRECTION VECTOR
    # ############
    
    # Relative wind speed in Airfoil's referential 'E'
    wE1 = wE[0]
    wE2 = wE[1]
    wE3 = wE[2]
    
    # Airfoil's transversal axis in carousel referential 'e'
#    eTe1 = e21
#    eTe2 = e22
#    eTe3 = e23
    
    # Lift axis ** Normed to we @>@> **
    eLe1 = - eTe2*we3 + eTe3*we2
    eLe2 = - eTe3*we1 + eTe1*we3
    eLe3 = - eTe1*we2 + eTe2*we1
       
    
    # AERODYNAMIC COEEFICIENTS
    # #################
    vT1 =          wE1
    vT2 = -lT*w3 + wE2
    vT3 =  lT*w2 + wE3
    
    
#    alpha = alpha0-wE3/wE1
#    beta = wE2/C.sqrt(wE1*wE1 + wE3*wE3)
    alpha = alpha0 + C.arctan2(-wE3,wE1)
    beta = C.arcsin(wE2/vKite)

    #NOTE: beta & alphaTail are compensated for the tail motion induced by
    #omega @>@>
#    alphaTail = alpha0-vT3/vT1
#    betaTail = vT2/C.sqrt(vT1*vT1 + vT3*vT3)
    alphaTail = alpha0 + C.arctan2(-vT3,vT1)
    betaTail = C.arcsin(vT2/vKite)
    
    dae['alpha(deg)'] = alpha*180/C.pi
    dae['alphaTail(deg)'] = alphaTail*180/C.pi
    dae['beta(deg)'] = beta*180/C.pi
    dae['betaTail(deg)'] = betaTail*180/C.pi

    # cL = cLA*alpha + cLe*elevator   + cL0
    # cD = cDA*alpha + cDA2*alpha*alpha + cDB2*beta*beta + cDe*elevator + cDr*aileron + cD0
    # cR = -rD*w1 + cRB*beta + cRAB*alphaTail*beta + cRr*aileron
    # cP = cPA*alphaTail + cPe*elevator  + cP0
    # cY = cYB*beta + cYAB*alphaTail*beta
    
    cL = cLA*alpha + cLe*elevator + cL0
    cD = cDA*alpha + cDA2*alpha*alpha + cDB2*beta*beta + cD0
    cR = -rD*w1 + cRB*betaTail + cRr*aileron + cRAB*alphaTail*betaTail
    cP = -pD*w2 + cPA*alphaTail + cPe*elevator + cP0
    cY = -yD*w3 + cYB*betaTail + cYAB*alphaTail*betaTail

    cD += 0.25*dae['r']*0.004/sref

    dae['cL'] = cL
    dae['cD'] = cD
    dae['L/D'] = cL/cD
    
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
     
    t1 =  0.5*rho*vKite2*bref*cR
    t2 =  0.5*rho*vKite2*cref*cP
    t3 =  0.5*rho*vKite2*bref*cY
    return (f1, f2, f3, t1, t2, t3)
