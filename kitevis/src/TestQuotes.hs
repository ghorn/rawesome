{-# OPTIONS_GHC -Wall #-}
-- {-# OPTIONS_GHC -ddump-splices #-}
{-# Language TemplateHaskell #-}

module TestQuotes where

import Control.Concurrent

import Quotes --( f, MyType(..) )

--import qualified Kite.CarouselState as CS
--import qualified Kite.Dcm as Dcm
--import qualified Kite.Xyz as KiteXyz

data Xyz = MkXyz { x_ :: Double
                 , y_ :: Double
                 , z_ :: Double
                 }
data Axyz = MkAxyz { a_ :: Double
                   , xyz_ :: Xyz
                   }

anAxyz :: Axyz
anAxyz = MkAxyz 7 (MkXyz 1 2 3)

increment :: Axyz -> Axyz
increment (MkAxyz a (MkXyz x y z)) = MkAxyz (a+1) (MkXyz (x + 0.1) (y + 0.2) (z + 0.3))

printVarInfo :: VarInfo -> IO ()
printVarInfo (VarInfo name type' mv) = do
  vals <- readMVar mv
  putStrLn $ "(" ++ show type' ++ ")  "++ name ++ ": " ++ show vals

go :: IO ()
go = do
  (receiveNewMessage, infos) <- $(setupTelem "position." ''Axyz)

  let printLog = mapM_ printVarInfo infos

      updateLoop 0 _ = return ()
      updateLoop n anAxyz' = do
        receiveNewMessage anAxyz'
        putStrLn ""
        printLog
        updateLoop (n-1::Int) (increment anAxyz')

  printLog
  updateLoop 4 anAxyz
  --woo <- $(f ''KiteXyz.Xyz)
  --woo <- $(f ''Dcm.Dcm)
  --woo <- $(f ''CS.CarouselState)
