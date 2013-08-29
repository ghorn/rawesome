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

module Main ( main ) where

import Data.Foldable ( toList )
import Data.Packed ( fromLists )
import Text.ProtocolBuffers.Basic ( uToString )

import SpatialMath

import MultiCarousel ( State(..), NiceKite(..), runMultiCarousel )
import qualified Carousel.Trajectory as PT
import qualified Carousel.Dae as PD
import qualified Carousel.DifferentialStates as PX

toNice :: PD.Dae -> NiceKite
toNice dae = NiceKite { nk_xyz = xyz
                      , nk_q'n'b = q'n'b
                      , nk_r'n0'a0 = r'n0'a0
                      , nk_r'n0't0 = r'n0't0
                      , nk_lineAlpha = 0.2 -- realToFrac $ fromMaybe 1 (CS.lineTransparency cs)
                      , nk_kiteAlpha = 1 -- realToFrac $ fromMaybe 1 (CS.kiteTransparency cs)
                      , nk_visSpan = 3 -- fromMaybe 1 (CS.visSpan cs)
                      }
  where
    daeX = PD.differentialStates dae
--    daeZ = PD.algebraicVars dae
--    daeU = PD.controls dae
--    daeP = PD.parameters dae

    x = PX.x daeX
    y = PX.y daeX
    z = PX.z daeX

    e11 = PX.e11 daeX
    e12 = PX.e12 daeX
    e13 = PX.e13 daeX

    e21 = PX.e21 daeX
    e22 = PX.e22 daeX
    e23 = PX.e23 daeX

    e31 = PX.e31 daeX
    e32 = PX.e32 daeX
    e33 = PX.e33 daeX

    delta = atan2 (PX.sin_delta daeX) (PX.cos_delta daeX)

    q'n'a = Quat (cos(0.5*delta)) 0 0 (sin(0.5*delta))

    q'a'b = quatOfDcm $ fromLists [ [e11, e12, e13]
                                  , [e21, e22, e23]
                                  , [e31, e32, e33]
                                  ]
    q'n'b = q'n'a * q'a'b

    rArm = Xyz 1.2 0 0 -- CS.rArm
    xyzArm = rArm + Xyz x y z
    xyz = rotVecByQuatB2A q'n'a xyzArm

    zt = 0 --CS.zt cs
    r'n0'a0 = rotVecByQuatB2A q'n'a rArm
    r'n0't0 = xyz + (rotVecByQuatB2A q'n'b $ Xyz 0 0 zt)

toState :: PT.Trajectory -> State
toState ptTraj = State nicekites messages
  where
    messages = map uToString $ toList $ PT.messages ptTraj
    nicekites = map toNice $ toList $ PT.traj ptTraj

main :: IO ()
main = runMultiCarousel "carousel" toState
