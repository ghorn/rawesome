{-# OPTIONS_GHC -Wall #-}

module Joy ( Joy(..)
           , setupJoy
           , getJoy
           ) where

--import Control.Concurrent
import qualified Graphics.UI.SDL as SDL
import qualified Graphics.UI.SDL.Joystick as JS
import Graphics.UI.SDL.Types
import GHC.Int ( Int16 )

data Joy = Joy { jsAxes :: [Double]
               , jsButtons :: [Bool]
               , jsHats :: [[Hat]]
               }
--main :: IO ()
--main = do
--  js <- setupJoy
--  update js 100
--  print "done"
------  SDL.quit

setupJoy :: IO Joystick
setupJoy = do
  SDL.init [SDL.InitJoystick]
  JS.open 0
  

normalize :: Int16 -> Double
normalize k
  | k < 0     = fromIntegral k / (abs (fromIntegral (minBound ::Int16)))
  | otherwise = fromIntegral k / fromIntegral (maxBound ::Int16)

getJoy :: Joystick -> IO Joy
getJoy js = do
  JS.update
  axes <- mapM (JS.getAxis js . fromIntegral) [0..(JS.axesAvailable js - 1)]
  buttons <- mapM (JS.getButton js . fromIntegral) [0..(JS.buttonsAvailable js - 1)]
  hats <- mapM (JS.getHat js . fromIntegral) [0..(JS.hatsAvailable js - 1)]
  return $ Joy { jsAxes = map normalize axes
               , jsButtons = buttons
               , jsHats = hats
               }

--update :: Joystick -> Int -> IO ()
--update _ 0 = return ()
--update js n = do
--  j <- getJoy js
----  putStrLn $ show n ++ "\t"++show (jsAxes j)
--  putStrLn $ show n ++ "\t"++show (jsHats j)
--  threadDelay 300000
--  update js (n-1)
