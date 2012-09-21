{-# OPTIONS_GHC -Wall #-}

module Model ( modelInteg
             ) where

import Prelude hiding (span)
import Data.Packed
import Data.Packed.ST
import Numeric.Container

forcesTorques :: Vector Double -> Vector Double -> (Double,Double,Double,Double,Double,Double)
forcesTorques state u = (f1,f2,f3,t1,t2,t3)
  where
    rho =    1.23      --  density of the air             --  [ kg/m^3]
    rA = 1.085 --(dixit Kurt)
    alpha0 = -0*pi/180 
    
    --TAIL LENGTH
    lT = 0.4
    
    --ROLL DAMPING
    rD = 1e-2 
    pD = 0*1e-3
    yD = 0*1e-3
    
    --WIND-TUNNEL PARAMETERS
    --Lift (report p. 67)
    cLA = 5.064
    
    cLe = -1.924
    
    cL0 = 0.239
    
    --Drag (report p. 70)
    cDA = -0.195
    cDA2 = 4.268
    cDB2 = 5
    {-
    cDe = 0.044
    cDr = 0.111
    -}
    cD0 = 0.026
    
    --Roll (report p. 72)
    cRB = -0.062
    cRAB = -0.271 
    cRr = -5.637e-1
    
    --Pitch (report p. 74)
    cPA = 0.293
    cPe = -4.9766e-1
    
    cP0 = 0.03
    
    --Yaw (report p. 76)
    cYB = 0.05
    cYAB = 0.229
    
    span = 0.96
    chord = 0.1
    
    -----------------------     model integ ---------------------------------------
    x =   state @> 0
    y =   state @> 1
    
    e11 = state @> 3
    e12 = state @> 4
    e13 = state @> 5
    
    e21 = state @> 6
    e22 = state @> 7
    e23 = state @> 8
    
    e31 = state @> 9
    e32 = state @> 10
    e33 = state @> 11
    
    dx =  state @> 12
    dy =  state @> 13
    dz =  state @> 14
    
    w1 =  state @> 15
    w2 =  state @> 16
    w3 =  state @> 17

    ddelta = state @> 19
    
    u1 = u @> 1
    u2 = u @> 2
    

    --------------- kinfile -------------
    dpE = fromList [ dx*e11 + dy*e12 + dz*e13 + ddelta*e12*rA + ddelta*e12*x - ddelta*e11*y
                   , dx*e21 + dy*e22 + dz*e23 + ddelta*e22*rA + ddelta*e22*x - ddelta*e21*y
                   , dx*e31 + dy*e32 + dz*e33 + ddelta*e32*rA + ddelta*e32*x - ddelta*e31*y
                   ]
    dp_carousel_frame = fromList [ dx - ddelta*y
                                 , dy + ddelta*rA + ddelta*x
                                 , dz
                                 ]
    
    ---------- more model_integ ----------------------
    -- EFFECTIVE WIND IN THE KITE`S SYSTEM :
    -- ---------------------------------------------------------------
    
    --Airfoil speed in carousel frame
    we1 = dp_carousel_frame @> 0
    we2 = dp_carousel_frame @> 1
    we3 = dp_carousel_frame @> 2
    
    vKite2 = dp_carousel_frame <.> dp_carousel_frame --Airfoil speed^2 
    vKite = sqrt(vKite2) --Airfoil speed
    
    -- CALCULATION OF THE FORCES :
    -- ---------------------------------------------------------------
    --
    --   FORCE ARE COMPUTED IN THE CAROUSEL FRAME !!!
    
    --Aero coeff.
    
            
    -- LIFT DIRECTION VECTOR
    -- -------------------------
    
    -- Relative wind speed in Airfoil's referential 'E'
    
    wE1 = dpE @> 0
    wE2 = dpE @> 1
    wE3 = dpE @> 2
    
    -- Airfoil's transversal axis in carousel referential 'e'
    eTe1 = e21
    eTe2 = e22
    eTe3 = e23
    
    
    -- Lift axis ** Normed to we !! **
    eLe1 = - eTe2*we3 + eTe3*we2
    eLe2 = - eTe3*we1 + eTe1*we3
    eLe3 = - eTe1*we2 + eTe2*we1
       
    
    -- AERODYNAMIC COEEFICIENTS
    -- ----------------------------------
    vT1 =          wE1
    vT2 = -lT*w3 + wE2
    vT3 =  lT*w2 + wE3
    
    
    alpha = alpha0-wE3/wE1
    
    --NOTE: beta & alphaTail are compensated for the tail motion induced by
    --omega !!
    betaTail = vT2/sqrt(vT1*vT1 + vT3*vT3)
    beta = wE2/sqrt(wE1*wE1 + wE3*wE3)
    alphaTail = alpha0-vT3/vT1
    
    -- cL = cLA*alpha + cLe*u2   + cL0
    -- cD = cDA*alpha + cDA2*alpha*alpha + cDB2*beta*beta + cDe*u2 + cDr*u1 + cD0
    -- cR = -rD*w1 + cRB*beta + cRAB*alphaTail*beta + cRr*u1
    -- cP = cPA*alphaTail + cPe*u2  + cP0
    -- cY = cYB*beta + cYAB*alphaTail*beta
    
    cL' = cLA*alpha + cLe*u2 + cL0
    cD' = cDA*alpha + cDA2*alpha*alpha + cDB2*beta*beta + cD0
    cR = -rD*w1 + cRB*betaTail + cRr*u1 + cRAB*alphaTail*betaTail
    cP = -pD*w2 + cPA*alphaTail + cPe*u2 + cP0
    cY = -yD*w3 + cYB*betaTail + cYAB*alphaTail*betaTail
    
    
    -- LIFT :
    -- ---------------------------------------------------------------
    cL = 0.2*cL'
    cD = 0.5*cD'
    fL1 =  rho*cL*eLe1*vKite/2.0
    fL2 =  rho*cL*eLe2*vKite/2.0
    fL3 =  rho*cL*eLe3*vKite/2.0
    
    -- DRAG :
    -- -----------------------------------------------------------
    fD1 = -rho*vKite*cD*we1/2.0
    fD2 = -rho*vKite*cD*we2/2.0 
    fD3 = -rho*vKite*cD*we3/2.0 
    
    
    
    -- FORCES (AERO)
    -- ---------------------------------------------------------------
    f1 = fL1 + fD1
    f2 = fL2 + fD2
    f3 = fL3 + fD3
    
    --f = f-fT
       
    -- TORQUES (AERO)
    -- ---------------------------------------------------------------
     
    t1 =  0.5*rho*vKite2*span*cR
    t2 =  0.5*rho*vKite2*chord*cP
    t3 =  0.5*rho*vKite2*span*cY

    
modelInteg :: Double -> Vector Double -> Vector Double -> (Vector Double, Vector Double)
modelInteg r state u = (sys, fromList [c, cdot, cddot])
  where
    --  PARAMETERS OF THE KITE :
    --  -----------------------------
    m =  0.626      --  mass of the kite               --  [ kg    ]
                 
    --   PHYSICAL CONSTANTS :
    --  -----------------------------
    g =    9.81      --  gravitational constant         --  [ m /s^2]
    
    --  PARAMETERS OF THE CABLE :
    --  -----------------------------
     
    --CAROUSEL ARM LENGTH
    rA = 1.085 --(dixit Kurt)

    zt = -0.01
    
    --INERTIA MATRIX (Kurt's direct measurements)
    j1 = 0.0163
    j31 = 0.0006
    j2 = 0.0078
    j3 = 0.0229
    
    --Carousel Friction & inertia
    jCarousel = 1e2
    cfric = 100

    
    -----------------------     model integ ---------------------------------------
    x =   state @> 0
    y =   state @> 1
    z =   state @> 2
    
    e11 = state @> 3
    e12 = state @> 4
    e13 = state @> 5
    
    e21 = state @> 6
    e22 = state @> 7
    e23 = state @> 8
    
    e31 = state @> 9
    e32 = state @> 10
    e33 = state @> 11
    
    dx =  state @> 12
    dy =  state @> 13
    dz =  state @> 14
    
    w1 =  state @> 15
    w2 =  state @> 16
    w3 =  state @> 17

    ddelta = state @> 19
    
    tc = u @> 0 --Carousel motor torque

    (f1,f2,f3,t1,t2,t3) = forcesTorques state u

    -- ATTITUDE DYNAMICS
    -- -----------------------------------------------------------
    
    -- DynFile          -- Call DAE
    mm :: Matrix Double
    mm = runSTMatrix $ do
      mm' <- newMatrix 0 8 8
      let writeMatrix' (row,col) = writeMatrix mm' (row-1) (col-1)
      writeMatrix' (1,1) $ jCarousel + m*rA*rA + m*x*x + m*y*y + 2*m*rA*x 
      writeMatrix' (1,2) $ -m*y 
      writeMatrix' (1,3) $ m*(rA + x) 
      writeMatrix' (1,4) $ 0 
      writeMatrix' (1,5) $ 0 
      writeMatrix' (1,6) $ 0 
      writeMatrix' (1,7) $ 0 
      writeMatrix' (1,8) $ 0
      
      writeMatrix' (2,1) $ -m*y 
      writeMatrix' (2,2) $ m 
      writeMatrix' (2,3) $ 0 
      writeMatrix' (2,4) $ 0 
      writeMatrix' (2,5) $ 0 
      writeMatrix' (2,6) $ 0 
      writeMatrix' (2,7) $ 0 
      writeMatrix' (2,8) $ x + zt*e31
      
      writeMatrix' (3,1) $ m*(rA + x) 
      writeMatrix' (3,2) $ 0 
      writeMatrix' (3,3) $ m 
      writeMatrix' (3,4) $ 0 
      writeMatrix' (3,5) $ 0 
      writeMatrix' (3,6) $ 0 
      writeMatrix' (3,7) $ 0 
      writeMatrix' (3,8) $ y + zt*e32
      
      writeMatrix' (4,1) $ 0 
      writeMatrix' (4,2) $ 0 
      writeMatrix' (4,3) $ 0 
      writeMatrix' (4,4) $ m 
      writeMatrix' (4,5) $ 0 
      writeMatrix' (4,6) $ 0 
      writeMatrix' (4,7) $ 0 
      writeMatrix' (4,8) $ z + zt*e33
      
      writeMatrix' (5,1) $ 0 
      writeMatrix' (5,2) $ 0 
      writeMatrix' (5,3) $ 0 
      writeMatrix' (5,4) $ 0 
      writeMatrix' (5,5) $ j1 
      writeMatrix' (5,6) $ 0 
      writeMatrix' (5,7) $ j31 
      writeMatrix' (5,8) $ -zt*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33)
      
      writeMatrix' (6,1) $ 0 
      writeMatrix' (6,2) $ 0 
      writeMatrix' (6,3) $ 0 
      writeMatrix' (6,4) $ 0 
      writeMatrix' (6,5) $ 0 
      writeMatrix' (6,6) $ j2 
      writeMatrix' (6,7) $ 0 
      writeMatrix' (6,8) $ zt*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33)
      
      writeMatrix' (7,1) $ 0 
      writeMatrix' (7,2) $ 0 
      writeMatrix' (7,3) $ 0 
      writeMatrix' (7,4) $ 0 
      writeMatrix' (7,5) $ j31 
      writeMatrix' (7,6) $ 0 
      writeMatrix' (7,7) $ j3 
      writeMatrix' (7,8) $ 0
      
      writeMatrix' (8,1) $ -zt*(e11*e23*x - e13*e21*x + e12*e23*y - e13*e22*y + zt*e11*e23*e31 - zt*e13*e21*e31 + zt*e12*e23*e32 - zt*e13*e22*e32) 
      writeMatrix' (8,2) $ x + zt*e31 
      writeMatrix' (8,3) $ y + zt*e32 
      writeMatrix' (8,4) $ z + zt*e33 
      writeMatrix' (8,5) $ -zt*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33) 
      writeMatrix' (8,6) $ zt*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) 
      writeMatrix' (8,7) $ 0 
      writeMatrix' (8,8) $ 0
      return mm'

    rhs :: Vector Double
    rhs = fromList
          [ tc - cfric*ddelta - f1*y + f2*(rA + x) + dy*m*(dx - 2*ddelta*y) - dx*m*(dy + 2*ddelta*rA + 2*ddelta*x) 
          , f1 + ddelta*m*(dy + ddelta*rA + ddelta*x) + ddelta*dy*m 
          , f2 - ddelta*m*(dx - ddelta*y) - ddelta*dx*m 
          , f3 - g*m 
          , t1 - w2*(j3*w3 + j31*w1) + j2*w2*w3 
          , t2 + w1*(j3*w3 + j31*w1) - w3*(j1*w1 + j31*w3) 
          , t3 + w2*(j1*w1 + j31*w3) - j2*w1*w2 
          , (w1 - ddelta*e13)*(e21*(zt*dx - zt2*e21*(w1 - ddelta*e13) + zt2*e11*(w2 - ddelta*e23)) + e22*(zt*dy - zt2*e22*(w1 - ddelta*e13) + zt2*e12*(w2 - ddelta*e23)) + zt*e33*(z*w1 + ddelta*e11*x + ddelta*e12*y + zt*e33*w1 + zt*ddelta*e11*e31 + zt*ddelta*e12*e32) + zt*e23*(dz + zt*e13*w2 - zt*e23*w1) + zt*e31*(w1 - ddelta*e13)*(x + zt*e31) + zt*e32*(w1 - ddelta*e13)*(y + zt*e32)) - dz*(dz + zt*e13*w2 - zt*e23*w1) - dx*(dx - zt*e21*(w1 - ddelta*e13) + zt*e11*(w2 - ddelta*e23)) - dy*(dy - zt*e22*(w1 - ddelta*e13) + zt*e12*(w2 - ddelta*e23)) - (zt*w1*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) + zt*w2*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33))*(w3 - ddelta*e33) - (w2 - ddelta*e23)*(e11*(zt*dx - zt2*e21*(w1 - ddelta*e13) + zt2*e11*(w2 - ddelta*e23)) + e12*(zt*dy - zt2*e22*(w1 - ddelta*e13) + zt2*e12*(w2 - ddelta*e23)) - zt*e33*(z*w2 + ddelta*e21*x + ddelta*e22*y + zt*e33*w2 + zt*ddelta*e21*e31 + zt*ddelta*e22*e32) + zt*e13*(dz + zt*e13*w2 - zt*e23*w1) - zt*e31*(w2 - ddelta*e23)*(x + zt*e31) - zt*e32*(w2 - ddelta*e23)*(y + zt*e32)) 
          ]
      where
        zt2 = zt*zt
    
    dRexp :: Matrix Double
    dRexp =
      fromLists [ [ e21*(w3 - ddelta*e33) - e31*(w2 - ddelta*e23) 
                  , e31*(w1 - ddelta*e13) - e11*(w3 - ddelta*e33) 
                  , e11*(w2 - ddelta*e23) - e21*(w1 - ddelta*e13) 
                  ]
                , [ e22*(w3 - ddelta*e33) - e32*(w2 - ddelta*e23) 
                  , e32*(w1 - ddelta*e13) - e12*(w3 - ddelta*e33) 
                  , e12*(w2 - ddelta*e23) - e22*(w1 - ddelta*e13) 
                  ]
                , [ e23*w3 - e33*w2 
                  , e33*w1 - e13*w3 
                  , e13*w2 - e23*w1
                  ]
                ]
    
    ddq_nu :: Vector Double
    ddq_nu = mm <\> rhs -- Backsolve DAE
    
    ----Extract dynamics
    --dddelta = ddq_nu(1)       -- Carousel acceleration
    --ddX = ddq_nu(2:4)         -- x,y,z acceleration
    --dw = ddq_nu(5:end-1)      -- Angular accelerations
    --nu = ddq_nu(end)          -- Algebraic state
    
{-    [dddelta, ddX, dw, nu] = takesV [1,3,3,1] ddq_nu -}
    [dddelta', ddX, dw, _] = takesV [1,3,3,1] ddq_nu
    
    --sys = [dx;dy;dz;reshape(dRexp,9,1);ddX;dw;ddelta;dddelta];
    sys = join [ fromList [dx,dy,dz]
               , flatten (trans dRexp)
               , ddX
               , dw
               , fromList [ddelta]
               , dddelta'
               ]

    (c,cdot,cddot) = (c',cdot', cddot')
      where
        dw1 = dw @> 0
        dw2 = dw @> 1
        {-
        dw3 = dw @> 2
        -}
        ddx = ddX @> 0
        ddy = ddX @> 1
        ddz = ddX @> 2
        dddelta = dddelta' @> 0
        
        c' =(x + zt*e31)**2/2 + (y + zt*e32)**2/2 + (z + zt*e33)**2/2 - r**2/2
        
        cdot' =dx*(x + zt*e31) + dy*(y + zt*e32) + dz*(z + zt*e33) + zt*(w2 - ddelta*e23)*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) - zt*(w1 - ddelta*e13)*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33)
        cddot' =(zt*w1*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) + zt*w2*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33))*(w3 - ddelta*e33) + dx*(dx + zt*e11*w2 - zt*e21*w1 - zt*ddelta*e11*e23 + zt*ddelta*e13*e21) + dy*(dy + zt*e12*w2 - zt*e22*w1 - zt*ddelta*e12*e23 + zt*ddelta*e13*e22) + dz*(dz + zt*e13*w2 - zt*e23*w1) + ddx*(x + zt*e31) + ddy*(y + zt*e32) + ddz*(z + zt*e33) - (w1 - ddelta*e13)*(e21*(zt*dx - zt**2*e21*(w1 - ddelta*e13) + zt**2*e11*(w2 - ddelta*e23)) + e22*(zt*dy - zt**2*e22*(w1 - ddelta*e13) + zt**2*e12*(w2 - ddelta*e23)) + zt*e33*(z*w1 + ddelta*e11*x + ddelta*e12*y + zt*e33*w1 + zt*ddelta*e11*e31 + zt*ddelta*e12*e32) + zt*e23*(dz + zt*e13*w2 - zt*e23*w1) + zt*e31*(w1 - ddelta*e13)*(x + zt*e31) + zt*e32*(w1 - ddelta*e13)*(y + zt*e32)) + (w2 - ddelta*e23)*(e11*(zt*dx - zt**2*e21*(w1 - ddelta*e13) + zt**2*e11*(w2 - ddelta*e23)) + e12*(zt*dy - zt**2*e22*(w1 - ddelta*e13) + zt**2*e12*(w2 - ddelta*e23)) - zt*e33*(z*w2 + ddelta*e21*x + ddelta*e22*y + zt*e33*w2 + zt*ddelta*e21*e31 + zt*ddelta*e22*e32) + zt*e13*(dz + zt*e13*w2 - zt*e23*w1) - zt*e31*(w2 - ddelta*e23)*(x + zt*e31) - zt*e32*(w2 - ddelta*e23)*(y + zt*e32)) + zt*(dw2 - dddelta*e23)*(e11*x + e12*y + e13*z + zt*e11*e31 + zt*e12*e32 + zt*e13*e33) - zt*(dw1 - dddelta*e13)*(e21*x + e22*y + e23*z + zt*e21*e31 + zt*e22*e32 + zt*e23*e33) - zt*dddelta*(e11*e23*x - e13*e21*x + e12*e23*y - e13*e22*y + zt*e11*e23*e31 - zt*e13*e21*e31 + zt*e12*e23*e32 - zt*e13*e22*e32)
