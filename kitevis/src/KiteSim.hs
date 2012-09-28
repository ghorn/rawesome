{-# OPTIONS_GHC -Wall #-}
{-# Language DoAndIfThenElse #-}

module Main where

import Data.Maybe ( fromMaybe )
import System.Random ( randomRs, mkStdGen)
import qualified System.ZMQ as ZMQ
import Control.Concurrent ( MVar, forkIO, modifyMVar_, newMVar, readMVar)
import Control.Monad ( forever )
import qualified Data.ByteString.Lazy as BL
import Data.Packed ( fromLists )
import Text.Printf ( printf )
import Text.ProtocolBuffers ( messageGet )

import qualified Kite.CarouselState as CS
import qualified Kite.Dcm as Dcm
import qualified Kite.Xyz as KiteXyz

import SpatialMath
import Vis
import Draw

data State = State { sTrails :: [[Xyz Double]]
                   , sCS :: Maybe CS.CarouselState
                   , sParticles :: [Xyz Double]
                   }

toNice :: CS.CarouselState -> (Xyz Double, Quat Double, Xyz Double, Xyz Double)
toNice cs = (xyz, q'n'b, r'n0'a0, r'n0't0)
  where
    x = KiteXyz.x $ CS.kiteXyz cs
    y = KiteXyz.y $ CS.kiteXyz cs
    z = KiteXyz.z $ CS.kiteXyz cs
    
    e11 = Dcm.r11 $ CS.kiteDcm cs
    e12 = Dcm.r12 $ CS.kiteDcm cs
    e13 = Dcm.r13 $ CS.kiteDcm cs

    e21 = Dcm.r21 $ CS.kiteDcm cs
    e22 = Dcm.r22 $ CS.kiteDcm cs
    e23 = Dcm.r23 $ CS.kiteDcm cs

    e31 = Dcm.r31 $ CS.kiteDcm cs
    e32 = Dcm.r32 $ CS.kiteDcm cs
    e33 = Dcm.r33 $ CS.kiteDcm cs

    delta = CS.delta cs

    q'nwu'ned = Quat 0 1 0 0

    q'n'a = Quat (cos(0.5*delta)) 0 0 (sin(-0.5*delta))

    q'aNWU'bNWU = quatOfDcmB2A $ fromLists [ [e11, e21, e31]
                                           , [e12, e22, e32]
                                           , [e13, e23, e33]
                                           ]
    q'a'b = q'nwu'ned * q'aNWU'bNWU * q'nwu'ned
    q'n'b = q'n'a * q'a'b
    q'n'aNWU = q'n'a * q'nwu'ned

    rArm = Xyz (CS.rArm cs) 0 0
    xyzArm = rArm + Xyz x y z
    xyz = rotVecByQuatB2A q'n'aNWU xyzArm

    zt = CS.zt cs
    r'n0'a0 = rotVecByQuatB2A q'n'a rArm
    r'n0't0 = xyz + (rotVecByQuatB2A q'n'b $ Xyz 0 0 (-zt))

drawFun :: State -> VisObject Double
----drawFun state = VisObjects $ [axes] ++ (map text [-5..5]) ++ [boxText, ac, plane,trailLines]
drawFun (State {sCS=Nothing}) = VisObjects []
drawFun state@(State {sCS=Just cs}) = VisObjects $ [axes, txt, ac, plane, trailLines, arm, line, points]
  where
    (pos@(Xyz px py pz), quat, r'n0'a0, r'n0't0) = toNice cs

    points = VisPoints (sParticles state) (Just 2) $ makeColor 1 1 1 0.5
    
    axes = VisAxes (0.5, 15) (Xyz 0 0 0) (Quat 1 0 0 0)
    arm  = VisLine [Xyz 0 0 0, r'n0'a0] $ makeColor 1 1 0 1
    line = VisLine [r'n0'a0, r'n0't0]   $ makeColor 0 1 1 1
    plane = VisPlane (Xyz 0 0 1) 1 (makeColor 1 1 1 1) (makeColor 0.2 0.3 0.32 1)
--    text k = Vis2dText "KITEVIS 4EVER" (100,500 - k*100*x) TimesRoman24 (makeColor 0 (0.5 + x'/2) (0.5 - x'/2) 1)
--      where
--        x' = realToFrac $ (x + 1)/0.4*k/5
--    boxText = Vis3dText "I'm a plane" (Xyz 0 0 (x-0.2)) TimesRoman24 (makeColor 1 0 0 1)
    ddelta = CS.ddelta cs

    (u1,u2,tc,wind_x) = (CS.u1 cs, CS.u2 cs, CS.tc cs, fromMaybe 0 (CS.wind_x cs))
    txt = VisObjects
          [ Vis2dText (printf "x: %.3f" px) (30,90) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "y: %.3f" py) (30,60) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "z: %.3f" pz) (30,30) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "RPM: %.3f" (ddelta*60/(2*pi))) (30,120) TimesRoman24 (makeColor 1 1 1 1)
--          , Vis2dText (printf "c:   %.3g" c    ) (30,150) TimesRoman24 (makeColor 1 1 1 1)
--          , Vis2dText (printf "c':  %.3g" cdot ) (30,180) TimesRoman24 (makeColor 1 1 1 1)
--          , Vis2dText (printf "c'': %.3g" cddot) (30,210) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "u1: %.3g \t(*180/pi = %.3f)" u1 (u1*180/pi)) (30,300) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "u2: %.3g \t(*180/pi = %.3f)" u2 (u2*180/pi)) (30,270) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "tc: %.3g" tc)                                (30,240) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "wind_x: %.3g" wind_x)                        (30,210) TimesRoman24 (makeColor 1 1 1 1)
          ]
    (ac,_) = drawAc pos quat

    trailLines = drawTrails (sTrails state)


particleBox :: Double
particleBox = 4

state0 :: State
state0 = State { sCS = Nothing
               , sTrails = [[],[],[]]
               , sParticles = take 100 $ randomRs (Xyz (-particleBox) (-particleBox) (-particleBox),
                                                   Xyz particleBox particleBox particleBox)
                              (mkStdGen 0)
               }

updateTrail :: [Xyz a] -> Xyz a -> [Xyz a]
updateTrail trail0 xyz
  | length trail0 < 65 = xyz:trail0
  | otherwise = take 65 (xyz:trail0)

boundParticle :: Xyz Double -> Xyz Double
boundParticle xyz@(Xyz x y z)
  | x >  particleBox = boundParticle (Xyz (x-2*particleBox) y z)
  | x < -particleBox = boundParticle (Xyz (x+2*particleBox) y z)
  | y >  particleBox = boundParticle (Xyz x (y-2*particleBox) z)
  | y < -particleBox = boundParticle (Xyz x (y+2*particleBox) z)
  | z >  particleBox = boundParticle (Xyz x y (z-2*particleBox))
  | z < -particleBox = boundParticle (Xyz x y (z+2*particleBox))
  | otherwise = xyz

updateState :: CS.CarouselState -> State -> IO State
updateState cs x0 =
  return $ State { sCS = Just cs
                 , sTrails = zipWith updateTrail trails0 trails
                 , sParticles = map (boundParticle . ((Xyz (ts*wind_x) 0 0) +)) (sParticles x0)
                 }
  where
    wind_x = fromMaybe 0 (CS.wind_x cs)
    trails0 = sTrails x0
    (pos,q,_,_) = toNice cs
    (_,trails) = drawAc pos q
    
sub :: MVar State -> IO ()
sub m = ZMQ.withContext 1 $ \context -> do
  ZMQ.withSocket context ZMQ.Sub $ \subscriber -> do
    ZMQ.connect subscriber "tcp://localhost:5563"
    ZMQ.subscribe subscriber "carousel"
    forever $ do
      addr <- ZMQ.receive subscriber []
      mre <- ZMQ.moreToReceive subscriber
      if mre
      then do
        msg <- ZMQ.receive subscriber []
        let cs = case messageGet (BL.fromChunks [msg]) of
              Left err -> error err
              Right (cs',_) -> cs'
        modifyMVar_ m (updateState cs)
      else return ()

ts :: Double
ts = 0.02

main :: IO ()
main = do
  m <- newMVar state0
  _ <- forkIO (sub m)
  
--  threadDelay 5000000
  let simFun _ _ = return ()
      df _ = fmap drawFun (readMVar m)
  simulateIO ts () df simFun
