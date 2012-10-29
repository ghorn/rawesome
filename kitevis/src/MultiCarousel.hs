{-# OPTIONS_GHC -Wall #-}
{-# Language DoAndIfThenElse #-}

module Main where

import Data.Foldable ( toList )
import Data.Maybe ( fromMaybe )
import qualified System.ZMQ as ZMQ
import Control.Concurrent ( MVar, forkIO, modifyMVar_, newMVar, readMVar)
import Control.Monad ( forever )
import qualified Data.ByteString.Lazy as BL
import Data.Packed ( fromLists )
import Text.ProtocolBuffers ( messageGet )
import Text.ProtocolBuffers.Basic ( uToString )

import qualified Kite.MultiCarousel as MC
import qualified Kite.CarouselState as CS
import qualified Kite.Dcm as Dcm
import qualified Kite.Xyz as KiteXyz

import SpatialMath
import Vis
import DrawAC
import ParseArgs ( getip )

type State = Maybe MC.MultiCarousel

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

    q'aNWU'bNWU = quatOfDcm $ fromLists [ [e11, e12, e13]
                                        , [e21, e22, e23]
                                        , [e31, e32, e33]
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

drawOneKite :: CS.CarouselState -> (VisObject Double, Double)
drawOneKite cs = (VisObjects [ac, arm, line], z)
  where
    (pos@(Xyz _ _ z), quat, r'n0'a0, r'n0't0) = toNice cs
    lineAlpha = realToFrac $ fromMaybe 1 (CS.lineTransparency cs)
    kiteAlpha = realToFrac $ fromMaybe 1 (CS.kiteTransparency cs)
    arm  = Line [Xyz 0 0 0, r'n0'a0] $ makeColor 1 1 0 lineAlpha
    line = Line [r'n0'a0, r'n0't0]   $ makeColor 0 1 1 lineAlpha

    (ac,_) = drawAc kiteAlpha pos quat

drawFun :: State -> VisObject Double
drawFun Nothing = VisObjects []
drawFun (Just ko) = VisObjects $ [axes,txt] ++ [plane] ++ kites
  where
    (kites,zs) = unzip $ map drawOneKite (toList (MC.css ko))
    z = maximum zs
--    points = Points (sParticles state) (Just 2) $ makeColor 1 1 1 0.5
    
    axes = Axes (0.5, 15)
--    arm  = Line [Xyz 0 0 0, r'n0'a0] $ makeColor 1 1 0 1
--    line = Line [r'n0'a0, r'n0't0]   $ makeColor 0 1 1 1
    plane = Trans (Xyz 0 0 planeZ') $ Plane (Xyz 0 0 1) (makeColor 1 1 1 1) (makeColor 0.2 0.3 0.32 planeAlpha)

    planeZ' = 2
    planeAlpha = 1*planeAlphaFade + (1 - planeAlphaFade)*0.2
    planeAlphaFade
      | z < planeZ' = 1
      | z < planeZ'+2 = realToFrac $ (planeZ'+2-z)/2
      | otherwise = 0
--    text k = 2dText "KITEVIS 4EVER" (100,500 - k*100*x) TimesRoman24 (makeColor 0 (0.5 + x'/2) (0.5 - x'/2) 1)
--      where
--        x' = realToFrac $ (x + 1)/0.4*k/5
--    boxText = 3dText "I'm a plane" (Xyz 0 0 (x-0.2)) TimesRoman24 (makeColor 1 0 0 1)
--    ddelta = CS.ddelta cs

--    (u1,u2,tc,wind_x) = (CS.u1 cs, CS.u2 cs, CS.tc cs, CS.wind_x cs)
    txt = VisObjects $
          zipWith (\s k -> Text2d (uToString s) (30,fromIntegral $ 30*k) TimesRoman24 (makeColor 1 1 1 1)) messages (reverse [1..length messages])
    messages = toList $ MC.messages ko


particleBox :: Double
particleBox = 4

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

updateState :: MC.MultiCarousel -> State -> IO State
updateState ko _ = return $ Just ko

sub :: String -> MVar State -> IO ()
sub ip m = ZMQ.withContext 1 $ \context -> do
  ZMQ.withSocket context ZMQ.Sub $ \subscriber -> do
    ZMQ.connect subscriber ip
    ZMQ.subscribe subscriber "multi-carousel"
    forever $ do
      _ <- ZMQ.receive subscriber []
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
  ip <- getip "multicarousel" "tcp://localhost:5563"

  m <- newMVar Nothing
  _ <- forkIO (sub ip m)
  
--  threadDelay 5000000
  let simFun _ _ = return ()
      df _ = fmap drawFun (readMVar m)
  simulateIO (Just ((1260,940),(1930,40))) "multi-carousel" ts () df simFun
