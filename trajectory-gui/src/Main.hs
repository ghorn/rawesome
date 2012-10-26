{-# OPTIONS_GHC -Wall #-}

--import Control.Concurrent
import Graphics.UI.Gtk

--import Data.Foldable ( toList )
--import Data.Maybe ( fromMaybe )
--import qualified System.ZMQ as ZMQ
--import Control.Concurrent ( MVar, forkIO, modifyMVar_, newMVar, readMVar)
--import Control.Monad ( forever )
import qualified Data.ByteString.Lazy as BL
--import Data.Packed ( fromLists )
import Text.ProtocolBuffers ( messagePut )
import Text.ProtocolBuffers.Basic ( uFromString )

import Control.Monad ( liftM )
import qualified Kite.MultiCarousel as MS
import qualified Kite.CarouselState as CS
import qualified Kite.Dcm as Dcm
import qualified Kite.Xyz as KiteXyz

import Data.Sequence ( fromList )
import System.ZMQ
import Data.ByteString.Char8 (pack)
--import Control.Concurrent (threadDelay)
--import Control.Monad (forever)

main :: IO ()
main = withContext 1 $ \context -> do
    withSocket context Pub $ \publisher -> do
      bind publisher "tcp://*:5563"
      main' publisher
--        forever $ do
--            send publisher (pack "A") [SndMore]
--            send publisher (pack "We don't want to see this") []
--            send publisher (pack "B") [SndMore]
--            send publisher (pack "We would like to see this") []
--            threadDelay $ 1 * 1000 * 1000

toCarousel :: [Double] -> CS.CarouselState
toCarousel x =
  CS.CarouselState
    { CS.kiteXyz = KiteXyz.Xyz (x !! 0) (x !! 1) (x !! 2)
    , CS.kiteDcm = Dcm.Dcm
                   (x !! 3) (x !! 4) (x !! 5)
                   (x !! 6) (x !! 7) (x !! 8)
                   (x !! 9) (x !! 10) (x !! 11)
    , CS.delta = x !! 18
    , CS.rArm = 1
    , CS.zt = 0
    , CS.messages = fromList []
    , CS.w0 = Nothing
    , CS.transparency = Nothing
    }

sendPlanes :: System.ZMQ.Socket a -> [([Double], t, t1)] -> (Int, Int) -> IO ()
sendPlanes publisher log' (k0,kf) = do
  let carousels = map (\(x,_,_) -> toCarousel x) $ take (kf - k0 + 1) $ drop k0 log'
      multiCarousel = MS.MultiCarousel (fromList carousels) (fromList [uFromString $ "range: "++show (k0,kf)])
  send publisher (pack "multi-carousel") [SndMore]
  send publisher (BL.toStrict $ messagePut multiCarousel) []
  return ()

main' :: System.ZMQ.Socket a -> IO ()
main' publisher = do
  logLines <- fmap lines $ readFile "testrun.txt"
  let log' :: [([Double],[Double],Double)]
      log' = map read logLines

      n = length log'
  
  _ <- initGUI
  mainWindow <- windowNew

  sc0 <- vScaleNewWithRange 0 (fromIntegral n - 1) 1
  sc1 <- vScaleNewWithRange 0 (fromIntegral n - 1) 1
  
  scaleSetValuePos sc0 PosTop
  scaleSetValuePos sc1 PosTop
  scaleSetDigits sc0 0
  scaleSetDigits sc1 0
  rangeSetUpdatePolicy sc0 UpdateDiscontinuous
  rangeSetUpdatePolicy sc1 UpdateDiscontinuous
  
  rangeSetValue sc0 0
  rangeSetValue sc1 (fromIntegral n-1)

  let update = do
        k0 <- liftM floor $ rangeGetValue sc0
        kf <- liftM floor $ rangeGetValue sc1
        sendPlanes publisher log' (k0,kf)
        
  _ <- afterRangeValueChanged sc0 update
  _ <- afterRangeValueChanged sc1 update
    
  button <- buttonNewWithLabel "save"
  _ <- afterClicked button $ do
    k0 <- liftM floor $ rangeGetValue sc0
    kf <- liftM floor $ rangeGetValue sc1
    let log'' = map (\(x,u,_) -> x++u) $ take (kf - k0 + 1) $ drop k0 log'
    putStrLn "writing file..."
    writeFile "out.txt" $ unlines $ map show log''
    putStrLn "finished writing file"

  layout <- do
      vb <- vBoxNew False 0
      hb <- hBoxNew False 0
      boxPackStart hb sc0 PackNatural 0
      boxPackStart hb sc1 PackGrow 0
      boxPackStart vb hb PackGrow 0
      boxPackStart vb button PackNatural 0
      return vb

  windowSetTitle mainWindow "S.A.R.A.H."
  windowSetDefaultSize mainWindow 400 400
  _ <- on mainWindow objectDestroy mainQuit
  containerAdd mainWindow layout
  widgetShowAll mainWindow

  mainGUI
