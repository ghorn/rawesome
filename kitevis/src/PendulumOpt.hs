{-# OPTIONS_GHC -Wall #-}
{-# Language DoAndIfThenElse #-}

module Main where

import Data.Foldable ( toList )
import qualified System.ZMQ as ZMQ
import Control.Concurrent ( MVar, forkIO, modifyMVar_, newMVar, readMVar)
import Control.Monad ( forever )
import qualified Data.ByteString.Lazy as BL
import Text.Printf ( printf )
import Text.ProtocolBuffers ( messageGet )

import qualified Kite.PendulumOpt as PO

import SpatialMath
import Vis
import DrawAC

data State = State { sTrails :: [Xyz Double]
                   , sEndTime :: Double
                   , sIter :: Int
                   }


drawFun :: Maybe State -> VisObject Double
----drawFun state = VisObjects $ [axes] ++ (map text [-5..5]) ++ [boxText, ac, plane,trailLines]
drawFun Nothing = VisObjects []
drawFun (Just state) = VisObjects $ [axes, txt, plane, trailLines]
  where
    axes = Trans (Xyz 0 (-0.1) 0) $ Axes (0.5, 15)
    plane = Trans (Xyz 0 0 1) $ Plane (Xyz 0 0 1) (makeColor 1 1 1 1) (makeColor 0.2 0.3 0.32 1)
    txt = VisObjects
          [ Text2d (printf "iteration: %d" (sIter state)) (30,60) TimesRoman24 (makeColor 1 1 1 1)
          , Text2d (printf "endTime: %.3f" (sEndTime state)) (30,30) TimesRoman24 (makeColor 1 1 1 1)
          ]
    trailLines = drawTrail (sTrails state) (\a -> makeColor a (1-a) 0 1)


updateState :: PO.PendulumOpt -> (Maybe State) -> IO (Maybe State)
updateState po _ =
  return $ Just $ State { sIter = fromIntegral $ PO.iters po
                        , sEndTime = PO.endTime po
                        , sTrails = zipWith3 Xyz (toList (PO.x po)) [0,0..] (toList (PO.z po))
                        }
    
sub :: MVar (Maybe State) -> IO ()
sub m = ZMQ.withContext 1 $ \context -> do
  ZMQ.withSocket context ZMQ.Sub $ \subscriber -> do
    ZMQ.connect subscriber "tcp://localhost:5563"
    ZMQ.subscribe subscriber "pendulum-opt"
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
  m <- newMVar Nothing
  _ <- forkIO (sub m)
  
  let simFun _ _ = return ()
      df _ = fmap drawFun (readMVar m)
  simulateIO (Just ((1260,940),(1930,40))) "pendulum optimization" ts () df simFun
