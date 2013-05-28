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

module DrawAC ( drawAc
              , drawTrails
              , drawTrail
              ) where

import SpatialMath
import Vis

drawAc :: Float -> Xyz Double -> Quat Double -> (VisObject Double, [Xyz Double])
drawAc alpha pos quat = (VisObjects $ wing ++ [htail,vtail,body], vtip:wingtips)
  where
--    axes = Trans pos $ RotQuat quat $ Axes (0.5, 15)
    spanW = 2.9
    arW = 9
    
    spanHRatio = 0.4
    arH = 4

    spanVRatio = 0.15
    arV = 2
    
    lengthToSpan = 0.7
    
    len = spanW*lengthToSpan
    width = 0.05
    
    wingPosOverLen = 0.2
    
    deltaWingTail = len*(1-wingPosOverLen)
    
    dcm = dcmOfQuat quat
    rotateTranslate = (pos +) . (rotVecByDcmB2A dcm)

    spanV = spanW*spanVRatio
    
    wingtips = map rotateTranslate [Xyz 0 (-spanW/2) 0, Xyz 0 (spanW/2) 0]
    vtip = rotateTranslate $ Xyz (-deltaWingTail) 0 (-spanV)
    
    wing = [ Quad
             (rotateTranslate $ Xyz ( chordW/2) ( spanW/2) 0)
             (rotateTranslate $ Xyz ( chordW/2) (-spanW/2) 0)
             (rotateTranslate $ Xyz (-chordW/2) (-spanW/2) 0)
             (rotateTranslate $ Xyz (-chordW/2) ( spanW/2) 0)
             $ makeColor 0 0 1 alpha
           , Quad
             (rotateTranslate $ Xyz (-chordW/2) ( spanW/2) (-0.01))
             (rotateTranslate $ Xyz (-chordW/2) (-spanW/2) (-0.01))
             (rotateTranslate $ Xyz ( chordW/2) (-spanW/2) (-0.01))
             (rotateTranslate $ Xyz ( chordW/2) ( spanW/2) (-0.01))
             $ makeColor 1 1 0 alpha
           ]
      where
        chordW = spanW/arW
    
        
    htail = Quad
            (rotateTranslate $ Xyz (-deltaWingTail + chordH/2) ( spanH/2) 0)
            (rotateTranslate $ Xyz (-deltaWingTail + chordH/2) (-spanH/2) 0)
            (rotateTranslate $ Xyz (-deltaWingTail - chordH/2) (-spanH/2) 0)
            (rotateTranslate $ Xyz (-deltaWingTail - chordH/2) ( spanH/2) 0)
            $ makeColor 0 0 1 alpha
      where
        spanH = spanW*spanHRatio
        chordH = spanH/arH
        
    vtail = Quad
            (rotateTranslate $ Xyz (-deltaWingTail + chordV/2) 0 (-spanV))
            (rotateTranslate $ Xyz (-deltaWingTail + chordV/2) 0 (     0))
            (rotateTranslate $ Xyz (-deltaWingTail - chordV/2) 0 (     0))
            (rotateTranslate $ Xyz (-deltaWingTail - chordV/2) 0 (-spanV))
            $ makeColor 1 1 0 alpha
      where
        chordV = spanV/arV
            
    body = Trans (rotateTranslate (Xyz (len/2-deltaWingTail) 0 0)) $ RotQuat quat $ 
           Ellipsoid (len/2, width/2, width/2) Solid (makeColor 0 0 1 alpha)

drawTrails :: [[Xyz a]] -> VisObject a
drawTrails xyzs = VisObjects $ zipWith drawTrail xyzs $ cycle [makeColor 0 0 1, makeColor 1 0 0, makeColor 0 1 0]

drawTrail :: [Xyz a] -> (Float -> Color) -> VisObject a
drawTrail trail mkCol = Line' $ zip trail (map mkCol (linspace 1 0 (length trail)))
  where
    linspace :: Fractional a => a -> a -> Int -> [a]
    linspace x0 xf n = map (\k -> x0 + (xf - x0) * (fromIntegral k) / (fromIntegral n-1)) [0..(n-1)]
