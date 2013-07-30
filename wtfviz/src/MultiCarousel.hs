-- Copyright 2012-2013 Greg Horn
--
-- This file is part of rawesome.
--
-- rawesome is free software: you can redistribute it and/or modify
-- it under the terms of the GNU Lesser General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- rawesome is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU Lesser General Public License for more details.
--
-- You should have received a copy of the GNU Lesser General Public License
-- along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

{-# OPTIONS_GHC -Wall #-}
{-# Language OverloadedStrings #-}
{-# Language CPP #-}

module MultiCarousel ( NiceKite(..)
                     , State(..)
                     , runMultiCarousel
                     ) where

#if OSX
import qualified System.ZMQ3 as ZMQ
#else
import qualified System.ZMQ as ZMQ
#endif
import qualified Control.Concurrent as CC
import Control.Monad ( forever, when )
import qualified Data.ByteString.Lazy as BL
import Text.ProtocolBuffers ( messageGet )
import Text.ProtocolBuffers ( ReflectDescriptor, Wire )

import SpatialMath
import qualified Xyz
import Vis
import DrawAC
import ParseArgs ( getip )

data State = State [NiceKite] [String]

data NiceKite = NiceKite { nk_xyz :: Xyz Double
                         , nk_q'n'b :: Quat Double
                         , nk_r'n0'a0 :: Xyz Double
                         , nk_r'n0't0 :: Xyz Double
                         , nk_lineAlpha :: Float
                         , nk_kiteAlpha :: Float
                         , nk_visSpan :: Double
                         }

drawOneKite :: Double -> NiceKite -> (VisObject Double, Double)
drawOneKite minLineLength niceKite
  | nk_lineAlpha niceKite > 0 = (VisObjects [ac, arm, line], z)
  | otherwise = (ac, z)
  where
    pos@(Xyz _ _ z) = nk_xyz niceKite
    quat = nk_q'n'b niceKite
    r'n0'a0 = nk_r'n0'a0 niceKite
    r'n0't0 = nk_r'n0't0 niceKite
    
    arm  = Line [Xyz 0 0 0, r'n0'a0] $ makeColor 1 1 0 (nk_lineAlpha niceKite)
    line = VisObjects [ line1 $ makeColor 1 0.2 0 (nk_lineAlpha niceKite)
                      , line2 $ makeColor 0 1 1 (nk_lineAlpha niceKite)
                      ]
      where
        line1 = Line [r'n0'a0, rMid] 
        line2 = Line [rMid, r'n0't0] 

        rMid = r'n0'a0 + fmap (* (1 - minLineLength/normDr)) dr
        dr = r'n0't0 - r'n0'a0
        normDr = Xyz.norm dr
        

    s = nk_visSpan niceKite
    ac = Trans pos $ Scale (s,s,s) ac'
    (ac',_) = drawAc (nk_kiteAlpha niceKite) (Xyz 0 0 0) quat

drawFun :: State -> VisObject Double
drawFun (State [] _) = VisObjects []
drawFun (State niceKites messages) = VisObjects $ [axes,txt] ++ [plane] ++ kites
  where
    minLineLength = minimum $ map lineLength niceKites
      where
        lineLength nk = Xyz.norm (nk_r'n0't0 nk - nk_r'n0'a0 nk)

    (kites,zs) = unzip $ map (drawOneKite minLineLength) niceKites
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
          zipWith (\s k -> Text2d s (30,fromIntegral $ 30*k) TimesRoman24 (makeColor 1 1 1 1)) messages (reverse [1..length messages])

withContext :: (ZMQ.Context -> IO a) -> IO a
#if OSX
withContext = ZMQ.withContext
#else
withContext = ZMQ.withContext 1
#endif

sub :: (ReflectDescriptor a, Wire a) => String -> (a -> State) -> String -> CC.MVar State -> IO ()
sub channel toState ip m = withContext $ \context -> do
#if OSX
  let receive = ZMQ.receive
#else
  let receive = flip ZMQ.receive []
#endif
  ZMQ.withSocket context ZMQ.Sub $ \subscriber -> do
    ZMQ.connect subscriber ip
    ZMQ.subscribe subscriber channel
    forever $ do
      _ <- receive subscriber
      mre <- ZMQ.moreToReceive subscriber
      when mre $ do
        msg <- receive subscriber
        case messageGet (BL.fromChunks [msg]) of
              Left err -> return ()
              Right (cs,_) -> CC.swapMVar m (toState cs) >> return ()

ts :: Double
ts = 0.02

runMultiCarousel :: (ReflectDescriptor a, Wire a) => String -> (a -> State) -> IO ()
runMultiCarousel channel toState = do
  (ip,followkite) <- getip "tcp://localhost:5563"
  when followkite $ putStrLn "wtfviz doesn't respect followkite flag, yo"
  putStrLn $ "using ip \""++ip++"\""

  m <- CC.newMVar (State [] [])
  _ <- CC.forkIO (sub channel toState ip m)
  
--  threadDelay 5000000
  let simFun _ _ = return ()
      df _ = fmap drawFun (CC.readMVar m)
  simulateIO (Just ((1260,940),(1930,40))) "multi-carousel" ts () df simFun
