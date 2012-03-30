{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

-- | RNAfold combinators, extracted for quickcheck

module BioInf.RNAfold.Combinators
  ( (~~#)
  , (#~~)
  ) where

import Data.Array.Repa.Index
import qualified Data.Vector.Fusion.Stream as S
import Data.List

import qualified ADP.Fusion as F
import qualified ADP.Fusion.Monadic as M
import qualified ADP.Fusion.Monadic.Internal as F
import ADP.Fusion.Monadic (makeLeft_MinRight)
import ADP.Fusion.Monadic.Internal (Box(..))



infixl 9 #~~, ~~#

-- | The structure on the left is a subword with size 2-28. The maximal size
-- could be 30 but since the two combinators are linked, 29,30 would fail
-- anyways.

(#~~) = makeLeft_MinRight (3,30) 1
{-# INLINE (#~~) #-}

-- | The structure on the right is a subword with size 2-30, however we inspect
-- the stack an reduce the maximal size.

(~~#) xs ys = Box mk step xs ys where
  minT = 8  -- minimal total size of region
  minC = 3  -- minimal number of nuc's on the right
  maxC = 32
  {-# INLINE mk #-}
  mk (z:.k:.j,a,b) = let (_:.i) = z
                         cnsmd = k-i -- consumed part
                         l = max k (j-maxC+cnsmd)
                     in return (z:.k:.l:.j,a,b)
  {-# INLINE step #-}
  step (z:.k:.l:.j,a,b)
    | l<=j-(max 0 $ minT - cnsmd) && l+minC<=j
    = return $ S.Yield (z:.k:.l:.j,a,b) (z:.k:.l+1:.j,a,b)
    | otherwise = return $ S.Done
    where cnsmd = k-i
          (_:.i) = z
{-# INLINE (~~#) #-}

