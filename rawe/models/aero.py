import casadi as C

def getWindAnglesFrom_v_bw_b(airspeed, v_bw_b):
    alpha =  C.arctan2( v_bw_b[2], v_bw_b[0] )
    beta  =  C.arcsin (-v_bw_b[1] / airspeed )
    return (alpha, beta)

def aeroForcesTorques(dae, conf, v_bw_n, v_bw_b, (w1,w2,w3), (eTe1, eTe2, eTe3), (aileron,elevator)):
    rho = conf['rho']
    alpha0 = conf['alpha0deg']*C.pi/180

    #ROLL DAMPING
    rD = conf['rD']
    pD = conf['pD']
    yD = conf['yD']

    #WIND-TUNNEL PARAMETERS
    #Lift (report p. 67)
    cLA = conf['cLA']

    cLe = conf['cLe']

    cL0 = conf['cL0']

    #Drag (report p. 70)
    cDA = conf['cDA']
    cDA2 = conf['cDA2']
    cDB2 = conf['cDB2']

    cD0 = conf['cD0']

    #Roll (report p. 72)
    cRB  = conf['cRB']
    cRAB = conf['cRAB']
    cRr  = conf['cRr']

    #Pitch (report p. 74)
    cPA = conf['cPA']
    cPe = conf['cPe']

    cP0 = conf['cP0']

    #Yaw (report p. 76)
    cYB = conf['cYB']
    cYAB = conf['cYAB']

    #TAIL LENGTH
    lT = conf['lT']

    sref = conf['sref']
    bref = conf['bref']
    cref = conf['cref']

    ##### more model_integ ###########
    # EFFECTIVE WIND IN THE KITE`S SYSTEM :
    # ###############################

    #Airfoil speed in carousel frame
    v_bw_n_x = v_bw_n[0]
    v_bw_n_y = v_bw_n[1]
    v_bw_n_z = v_bw_n[2]

    vKite2 = C.mul(v_bw_n.trans(), v_bw_n) #Airfoil speed^2
    vKite = C.sqrt(vKite2) #Airfoil speed
    dae['airspeed'] = vKite

    # LIFT DIRECTION VECTOR
    # ############

    # Relative wind speed in Airfoil's referential 'E'
    v_bw_b_x = v_bw_b[0]
    v_bw_b_y = v_bw_b[1]
    v_bw_b_z = v_bw_b[2]

    # Airfoil's transversal axis in carousel referential 'e'
#    eTe1 = e21
#    eTe2 = e22
#    eTe3 = e23

    # Lift axis ** Normed to we @>@> **
    eLe1 = - eTe2*v_bw_n_z + eTe3*v_bw_n_y
    eLe2 = - eTe3*v_bw_n_x + eTe1*v_bw_n_z
    eLe3 = - eTe1*v_bw_n_y + eTe2*v_bw_n_x


    # AERODYNAMIC COEEFICIENTS
    # #################
    v_tailw_b = C.veccat([          v_bw_b_x,
                           -lT*w3 + v_bw_b_y,
                            lT*w2 + v_bw_b_z])
    v_tailw_b_x = v_tailw_b[0]
    v_tailw_b_y = v_tailw_b[1]
    v_tailw_b_z = v_tailw_b[2]


    #NOTE: beta & alphaTail are compensated for the tail motion induced by
    #omega @>@>
    if conf['alpha_beta_computation'] == 'first_order':
        alpha = alpha0-v_bw_b_z/v_bw_b_x
        beta = v_bw_b_y/v_bw_b_x
#        beta = v_bw_b_y/C.sqrt(v_bw_b_x*v_bw_b_x + v_bw_b_z*v_bw_b_z)
        alphaTail = alpha0-v_tailw_b_z/v_tailw_b_x
        betaTail = v_tailw_b_y/v_tailw_b_x
#        betaTail = v_tailw_b_y/C.sqrt(v_tailw_b_x*v_tailw_b_x + v_tailw_b_z*v_tailw_b_z)

    elif conf['alpha_beta_computation'] == 'closed_form':
        (alpha, beta)          = getWindAnglesFrom_v_bw_b(vKite, v_bw_b)
        (alphaTail, betaTail)  = getWindAnglesFrom_v_bw_b(vKite, v_tailw_b)
        alpha     += alpha0
        alphaTail += alpha0
    else:
        raise ValueError('config "alpha_beta_compuation" value '+str(conf['alpha_beta_computation'])+' not recognized, use "first_order" or "closed_form"')

    dae['alpha_deg'] = alpha*180/C.pi
    dae['alphaTail_deg'] = alphaTail*180/C.pi
    dae['beta_deg'] = beta*180/C.pi
    dae['betaTail_deg'] = betaTail*180/C.pi

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

    dae['cL'] = cL
    dae['cD'] = cD
    dae['cD_tether'] = 0.25*dae['r']*0.001/sref
    dae['L_over_D'] = cL/cD
    cD = dae['cD'] + dae['cD_tether']
    dae['L_over_D_with_tether'] = cL/cD

    # LIFT :
    # ###############################
    dae['fL'] = sref*rho*cL*vKite2/2.0
    fL1 =  sref*rho*cL*eLe1*vKite/2.0
    fL2 =  sref*rho*cL*eLe2*vKite/2.0
    fL3 =  sref*rho*cL*eLe3*vKite/2.0

    # DRAG :
    # #############################
    dae['fD'] = sref*rho*vKite*cD*2.0
    fD1 = -sref*rho*vKite*cD*v_bw_n_x/2.0
    fD2 = -sref*rho*vKite*cD*v_bw_n_y/2.0
    fD3 = -sref*rho*vKite*cD*v_bw_n_z/2.0

    # FORCES (AERO)
    # ###############################
    f1 = fL1 + fD1
    f2 = fL2 + fD2
    f3 = fL3 + fD3

    #f = f-fT

    # TORQUES (AERO)
    # ###############################

    t1 =  0.5*rho*vKite2*sref*bref*cR
    t2 =  0.5*rho*vKite2*sref*cref*cP
    t3 =  0.5*rho*vKite2*sref*bref*cY

    dae['aero_fx'] = f1
    dae['aero_fy'] = f2
    dae['aero_fz'] = f3
    dae['aero_mx'] = t1
    dae['aero_my'] = t2
    dae['aero_mz'] = t3
    return (f1, f2, f3, t1, t2, t3)

#def tetherDragLumped( conf, r_n2b_n, v_bn_n, get_v_wn_n, tetherLength, N=10 ):
#    Cd = 1.1 # meh
#
#    forceTotal = C.DMatrix([0.0, 0.0, 0.0])
#    sumMe = 0
#    for k in range(N):
#        scaling = (k + 0.5)/N
#        sumM += scaling
#        r_n2t_n = r_n2b_n*scaling # local tether position
#        v_tn_n  =  v_bn_n*scaling
#        v_tw_b = v_tn_n - get_v_wn_n(r_n2t_n) # local tether velocity in wind frame
#        localV2 = C.mul(v_tw_b.T, v_tw_b)
#        localV = C.sqrt(localV2)
#
#        deltaH = scaling*r_n2b_n
#        hProjected = C.sqrt( C.mul(deltaH.T, deltaH) - C.mul(deltaH.T, v_tw_b)/localV
#
#        localForce = 0.5*rho*localV2*Cd*hProjected
#    raise Exception('summe: ',sumMe)
