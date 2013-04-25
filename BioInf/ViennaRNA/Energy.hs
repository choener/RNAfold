{-# LANGUAGE PatternGuards #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE BangPatterns #-}

module BioInf.ViennaRNA.Energy where

import Data.Vector.Fusion.Stream.Monadic as SM
import qualified Data.Vector.Unboxed as VU
import Control.Lens
import Data.Array.Repa.Index
import Prelude as P
import qualified Data.Map as M

import Data.PrimitiveArray as PA hiding ((!))
import Data.PrimitiveArray.Zero as PA
import qualified Data.PrimitiveArray as PA
import Biobase.Turner
import Biobase.Vienna
import Biobase.Primary

import BioInf.ViennaRNA.Signature

import Debug.Trace



mfe :: Monad m => Signature m Deka Deka
mfe = (hairpin,interior,multi,blockStem,blockUnpair,compsBR,compsBC,structW,structCS,structWS,structOpen,h) where
  hairpin ener l lp xs rp r
      | len <= 6
      , Just e <- (l `VU.cons` xs `VU.snoc` r) `M.lookup` _hairpinLookup ener = e
      | len <   3 = huge
      | len ==  3 = (ener^.hairpinL) VU.! len + tAU
      | len < 31  = (ener^.hairpinL) VU.! len + ener^.hairpinMM!(Z:.l:.r:.lp:.rp)
      | otherwise = huge
      where
        !len = VU.length xs
        !tAU  = if (l,r) == (nC,nG) || (l,r) == (nG,nC) then Deka 0 else ener^.termAU
  interior ener l ls li w ri rs r
      | lls==0 && lrs==0  -- stack
      = w + _stack ener ! (Z:.l:.r:.ri:.li) -- left, right, right inner, left inner
      | lls==1 && lrs==0 || lls==0 && lrs==1  -- stack with slip
      = w + _stack ener ! (Z:.l:.r:.ri:.li) + _bulgeL ener VU.! 1
      | lls==1 && lrs==1
      = w + _iloop1x1 ener ! (Z:.l:.r:.ri:.li:.lH:.rL)
      | lls==1 && lrs==2
      = w + _iloop2x1 ener ! (Z:.l:.r:.ri:.li:.lH:.rH:.rL)
      | lls==2 && lrs==1
      = w + _iloop2x1 ener ! (Z:.l:.r:.ri:.li:.rH:.lH:.lL)
      | lls==2 && lrs==2
      = w + _iloop2x2 ener ! (Z:.l:.r:.ri:.li:.lH:.lL:.rH:.rL)
      | min lls lrs == 2 && max lls lrs == 3
      = w + _iloop2x3MM ener ! (Z:.l:.r:.lH:.lL) + _iloop2x3MM ener ! (Z:.ri:.li:.rL:.rH) + _iloopL ener VU.! 5 + _ninio ener
      | lls==0 && lrs > 1 && lrs <= 30
      = w + tAU + _bulgeL ener VU.! lrs + tUA
      | lrs==0 && lls > 1 && lls <= 30
      = w + tAU + _bulgeL ener VU.! lls + tUA
      | lrs==1 && lls > 2 && lls <= 30
      = w + _iloop1xnMM ener ! (Z:.li:.ri:.lL:.rH) + _iloop1xnMM ener ! (Z:.r:.l:.rL:.lH) + _iloopL ener VU.! lls + min (_ninio ener *. (lls-1)) (_maxNinio ener)
      | lls==1 && lrs > 2 && lrs <= 30
      = w + _iloop1xnMM ener ! (Z:.li:.ri:.lL:.rH) + _iloop1xnMM ener ! (Z:.r:.l:.rL:.lH) + _iloopL ener VU.! lrs + min (_ninio ener *. (lrs-1)) (_maxNinio ener)
      | lls>0 && lrs>0 && lls+lrs <= 30 -- TODO missing support for length constraints ?
      = w + _iloopMM ener ! (Z:.l:.r:.lH:.rL) + _iloopMM ener ! (Z:.ri:.li:.rH:.lL) + _iloopL ener VU.! (lls+lrs) + min (_ninio ener *. (abs $ lls - lrs)) (_maxNinio ener)
      | otherwise = huge -- NOTE later on, we should never get this score
      where
        !lls = VU.length ls
        !lrs = VU.length rs
        !tAU = if (l,r)   `P.elem` [(nC,nG), (nG,nC)] then Deka 0 else ener^.termAU
        !tUA = if (li,ri) `P.elem` [(nC,nG), (nG,nC)] then Deka 0 else ener^.termAU
        lH = VU.unsafeHead ls
        lL = VU.unsafeLast ls
        rH = VU.unsafeHead rs
        rL = VU.unsafeLast rs
  multi ener l li b c ri r
    = b + c + _multiMM ener ! (Z:.r:.l:.ri:.li) + _multiHelix ener + _multiOffset ener where
  blockStem ener lo l s r ro
    = s + _multiMM ener ! (Z:.l:.r:.lo:.ro) + _multiHelix ener
  blockUnpair ener c b
    = b + _multiNuc ener
  compsBR ener b reg
    = let Deka nuc = _multiNuc ener in b + (Deka $ nuc * (VU.length reg))
  compsBC ener b c
    = b + c
  structW ener w
    = w
  structCS ener c w
    = w
  structWS ener w s
    = w + s
  structOpen ener r
    = 0
  h = foldl' min huge
  {-# INLINE hairpin #-}
  {-# INLINE interior #-}
  {-# INLINE multi #-}
  {-# INLINE blockStem #-}
  {-# INLINE blockUnpair #-}
  {-# INLINE compsBR #-}
  {-# INLINE compsBC #-}
  {-# INLINE structW #-}
  {-# INLINE structCS #-}
  {-# INLINE structWS #-}
  {-# INLINE structOpen #-}
  {-# INLINE h #-}
{-# INLINE mfe #-}

huge = Deka 999999
{-# INLINE huge #-}

infixl 8 !
(!) = (PA.!)
{-# INLINE (!) #-}

(*.) :: Deka -> Int -> Deka
(Deka k) *. n = Deka $ k*n
{-# INLINE (*.) #-}

