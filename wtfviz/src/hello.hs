{-# OPTIONS_GHC -Wall #-}
{-# Language TemplateHaskell #-}

module Main where

import Control.Monad ( zipWithM_ )
import qualified Data.Sequence as Seq
import qualified Data.Foldable as F
import Control.Concurrent ( forkIO, readMVar, threadDelay )
import Graphics.UI.Gtk ( AttrOp( (:=) ) )
import qualified Graphics.UI.Gtk as G
import qualified Graphics.UI.Gtk as Gtk
import qualified Graphics.UI.Gtk.ModelView as New
import qualified Graphics.UI.Gtk.OpenGL as GtkGL

import System.Glib.Signals (on)
import Data.List ( isPrefixOf )
import Data.Char ( toLower )

import Graphics.Rendering.OpenGL as GL
 
import Quotes --( f, MyType(..) )

data Xyz = MkXyz { x_ :: Double
                 , y_ :: Double
                 , z_ :: Double
                 }
data Axyz = MkAxyz { a_ :: Double
                   , xyz_ :: Xyz
                   }

data VarInfo' = VarInfo' { viName :: String
                         , viPc :: PContainer
                         , viMarked :: Bool
                         }

anAxyz :: Axyz
anAxyz = MkAxyz 7 (MkXyz 1 2 3)

increment :: Axyz -> Axyz
increment (MkAxyz a (MkXyz x y z)) = MkAxyz (a+1) (MkXyz (x + 0.1) (y + 0.2) (z + 0.3))

updateLoop :: Int -> Axyz -> (Axyz -> IO ()) -> IO ()
updateLoop n anAxyz' receiveNewMessage = do
  receiveNewMessage anAxyz'
  putStrLn $ "update " ++ show n
  threadDelay 1000000
  updateLoop (n+1::Int) (increment anAxyz') receiveNewMessage

main :: IO ()
main = do
  (receiveNewMessage, infos) <- $(setupTelem "position" ''Axyz)

  _ <- forkIO $ updateLoop 0 anAxyz receiveNewMessage

  _ <- G.initGUI

  win <- G.windowNew
  _ <- G.onDestroy win G.mainQuit

  -- create a new tree model
  model <- G.listStoreNew $ map (\(VarInfo st pc) -> VarInfo' st pc False) infos
  view <- New.treeViewNewWithModel model

  New.treeViewSetHeadersVisible view True

  -- add three columns
  col1 <- New.treeViewColumnNew
  col2 <- New.treeViewColumnNew
  col3 <- New.treeViewColumnNew

  New.treeViewColumnSetTitle col1 "name"
  New.treeViewColumnSetTitle col2 "Int column"
  New.treeViewColumnSetTitle col3 "show?"

  renderer1 <- New.cellRendererTextNew
  renderer2 <- New.cellRendererTextNew
  renderer3 <- New.cellRendererToggleNew

  New.cellLayoutPackStart col1 renderer1 True
  New.cellLayoutPackStart col2 renderer2 True
  New.cellLayoutPackStart col3 renderer3 True

  New.cellLayoutSetAttributes col1 renderer1 model $ \row -> [ New.cellText := viName row ]
  New.cellLayoutSetAttributes col2 renderer2 model $ \_ -> [ New.cellText := show "waaa" ]
  New.cellLayoutSetAttributes col3 renderer3 model $ \row -> [ New.cellToggleActive := viMarked row ]

  _ <- New.treeViewAppendColumn view col1
  _ <- New.treeViewAppendColumn view col2
  _ <- New.treeViewAppendColumn view col3

  -- update the model when the toggle buttons are activated
  _ <- on renderer3 G.cellToggled $ \pathStr -> do
    let (i:_) = G.stringToTreePath pathStr
    val <- G.listStoreGetValue model i
    G.listStoreSetValue model i val { viMarked = not (viMarked val) }


  -- enable interactive search
  New.treeViewSetEnableSearch view True
  New.treeViewSetSearchEqualFunc view $ Just $ \str iter -> do
    (i:_) <- G.treeModelGetPath model iter
    row <- G.listStoreGetValue model i
    return (map toLower str `isPrefixOf` map toLower (viName row))

  G.containerAdd win view
  G.widgetShowAll win

  setupGlstuff infos
  
  G.mainGUI


setupGlstuff :: [VarInfo] -> IO ()
setupGlstuff infos = do
--  _ <- Gtk.initGUI
 
  -- Initialise the Gtk+ OpenGL extension
  -- (including reading various command line parameters)
  _ <- GtkGL.initGL
 
  -- We need a OpenGL frame buffer configuration to be able to create other
  -- OpenGL objects.
  glconfig <- GtkGL.glConfigNew [GtkGL.GLModeRGBA,
                                 GtkGL.GLModeDepth,
                                 GtkGL.GLModeDouble]
 
  -- Create an OpenGL drawing area widget
  glCanvas <- GtkGL.glDrawingAreaNew glconfig
 
  _ <- Gtk.widgetSetSizeRequest glCanvas 250 250
 
  -- Initialise some GL setting just before the canvas first gets shown
  -- (We can't initialise these things earlier since the GL resources that
  -- we are using wouldn't heve been setup yet)
  _ <- Gtk.onRealize glCanvas $ GtkGL.withGLDrawingArea glCanvas $ \_ -> do
    clearColor $= (Color4 0.0 0.0 0.0 0.0)
    matrixMode $= Projection
    loadIdentity
    ortho 0.0 1.0 0.0 1.0 (-1.0) 1.0
    depthFunc $= Just Less
    drawBuffer $= BackBuffers
 
  -- Set the repaint handler
  _ <- Gtk.onExpose glCanvas $ \_ -> do
    GtkGL.withGLDrawingArea glCanvas $ \glwindow -> do
      GL.clear [GL.DepthBuffer, GL.ColorBuffer]
      display infos
      GtkGL.glDrawableSwapBuffers glwindow
    return True
 
  -- Setup the animation
  _ <- Gtk.timeoutAddFull (do
      Gtk.widgetQueueDraw glCanvas
      return True)
    Gtk.priorityDefaultIdle animationWaitTime
 
  --------------------------------
  -- Setup the rest of the GUI:
  --
  window <- Gtk.windowNew
  _ <- Gtk.onDestroy window Gtk.mainQuit
  _ <- Gtk.set window [ Gtk.containerBorderWidth := 8,
                        Gtk.windowTitle := "Gtk2Hs + HOpenGL demo" ]
 
  vbox <- Gtk.vBoxNew False 4
  _ <- Gtk.set window [ Gtk.containerChild := vbox ]
 
  label <- Gtk.labelNew (Just "Gtk2Hs using OpenGL via HOpenGL!")
  button <- Gtk.buttonNewWithLabel "Close"
  _ <- Gtk.onClicked button Gtk.mainQuit
  Gtk.set vbox [ Gtk.containerChild := glCanvas,
                 Gtk.containerChild := label,
                 Gtk.containerChild := button ]
 
  Gtk.widgetShowAll window
--  Gtk.mainGUI
 
--data PContainer = PCDouble (MVar (ContainerType Double))
--                | PCFloat (MVar (ContainerType Float))
--                | PCInt32 (MVar (ContainerType P'.Int32))
--                | PCInt64 (MVar (ContainerType P'.Int64))
--                | PCWord32 (MVar (ContainerType P'.Word32))
--                | PCWord64 (MVar (ContainerType P'.Word64))
--                | PCBool (MVar (ContainerType Bool))
--                | PCUtf8 (MVar (ContainerType P'.Utf8))
--                | PCByteString (MVar (ContainerType (BSL.ByteString)))
--
--data VarInfo = VarInfo String PContainer

linspace :: Fractional b => b -> b -> Int -> [b]
linspace x0 x1 n =
  map ((+ x0) . (((x1-x0)/((realToFrac n)-1)) *) . realToFrac) [0..(n-1)]

drawLine :: (Fractional a, Real a) => Seq.Seq a -> IO ()
drawLine seq' = do
  let maxval = F.foldl' (\acc x -> max acc (abs x)) (1e-12) seq'
      num = Seq.length seq'
      xs = linspace (-1) 1 num
      ys = map ((/ (realToFrac maxval)) . realToFrac) $ F.toList seq'

--  putStrLn $ "\nxs: " ++ show xs
--  putStrLn $ "ys: " ++ show ys
  renderPrimitive LineStrip $ do
    zipWithM_ (\x y -> vertex (Vertex3 x y 0.0 :: Vertex3 GLfloat)) xs ys

-- Draw the OpenGL polygon.
display :: [VarInfo] -> IO ()
display infos = do
--  let printLog = mapM_ printVarInfo infos
--  printLog
  loadIdentity
  color (Color3 1 1 1 :: Color3 GLfloat)
  -- Instead of glBegin ... glEnd there is renderPrimitive.
  let drawOne (VarInfo name (PCDouble mv)) = do
        s <- readMVar mv
        putStrLn $ "trying to draw " ++ name
        drawLine s
  mapM_ drawOne infos
--  renderPrimitive Polygon $ do
--    vertex (Vertex3 0.25 0.25 0.0 :: Vertex3 GLfloat)
--    vertex (Vertex3 0.75 0.25 0.0 :: Vertex3 GLfloat)
--    vertex (Vertex3 0.75 0.75 0.0 :: Vertex3 GLfloat)
--    vertex (Vertex3 0.25 0.75 0.0 :: Vertex3 GLfloat)
 
animationWaitTime :: Int
animationWaitTime = 3
