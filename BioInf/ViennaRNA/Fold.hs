{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE BangPatterns #-}

module BioInf.ViennaRNA.Fold where

import Data.Vector.Fusion.Util (Id (..))
import Data.Vector.Fusion.Stream.Monadic as SM
import qualified Data.Vector.Unboxed as VU
import Data.Array.Repa.Index
import Control.Monad
import Control.Monad.ST
import System.IO.Unsafe
import Prelude as P hiding (Maybe(..))
import Data.Strict.Maybe
import Data.Strict.Tuple

import Biobase.Secondary.Diagrams
import Data.Array.Repa.Index.Subword
import ADP.Fusion
import ADP.Fusion.Table
import Biobase.Vienna
import Biobase.Primary
import Data.PrimitiveArray as PA hiding ((!))
import Data.PrimitiveArray.Zero as PA

import BioInf.ViennaRNA.Signature
import BioInf.ViennaRNA.Energy



basepairing :: Primary -> Subword -> Bool
basepairing inp (Subword(i:.j)) = i+1<j && f (inp VU.! i) (inp VU.! (j-1)) where
  f l r =  l==nC && r==nG
        || l==nG && r==nC
        || l==nA && r==nU
        || l==nU && r==nA
        || l==nG && r==nU
        || l==nU && r==nG
  {-# INLINE f #-}
{-# INLINE basepairing #-}

structureConstrains :: Maybe D1Secondary -> Subword -> Bool
structureConstrains Nothing         !_               = True
structureConstrains !(Just (D1S c)) (Subword (i:.j)) = (i<j) && (VU.unsafeIndex c i == j-1)
{-# INLINE structureConstrains #-}

structC :: Primary -> Subword -> Bool
structC inp (Subword(i:.j)) = VU.length inp == j
{-# INLINE structC #-}

-- TODO need to fix sized regions, then we are good to go -- performance-wise
--
-- TODO backtracking
--
-- TODO struct table
--
-- TODO restrict structs to a linear band

gRNAfold ener (hairpin,interior,multi,blockStem,blockUnpair,compsBR,compsBC,structW,structCS,structWS,structOpen,h) weak block comps struct cs inp =
  ( weak ,
    hairpin  ener <<< c % pr % hr % pl % c             |||
    interior ener <<< c % ir % pr % weak % pl % ir % c |||
    multi    ener <<< c % pl % block % comps % pl % c `check` (basepairing inp) `check` (structureConstrains cs) ... h
  , block ,
    blockStem   ener <<< pl % c % weak % c % pr |||
    blockUnpair ener <<< c % block              ... h
  , comps ,
    compsBR ener <<< block % r     |||
    compsBC ener <<< block % comps ... h
  , struct ,
--    structW  ener <<< weak          |||       -- TODO peak left/right with default ; not needed anymore
    structCS ener <<< c % struct    |||
    structWS ener <<< weak % struct |||       -- peak here for weak, too
    structOpen ener <<< r           `check` (structC inp) ... h
  ) where c = chr inp
          r = region inp
          pr = peekR inp
          pl = peekL  inp
          hr = sregion 3 30 inp
          ir = sregion 0 20 inp
          {-# INLINE c #-}
          {-# INLINE r #-}
          {-# INLINE pr #-}
          {-# INLINE pl #-}
          {-# INLINE hr #-}
          {-# INLINE ir #-}
{-# INLINE gRNAfold #-}



pretty :: Monad m => Signature m String (SM.Stream m String)
pretty = (hairpin,interior,multi,blockStem,blockUnpair,compsBR,compsBC,structW,structCS,structWS,structOpen,h) where
  hairpin     _ _ _ r _ _ = "(" P.++ (P.replicate (VU.length r) '.') P.++ ")"
  interior    _ _ l _ w _ r _ = "(" P.++ (P.replicate (VU.length l) '.') P.++ w P.++ (P.replicate (VU.length r) '.') P.++ ")"
  multi       _ _ _ b c _ _ = "(" P.++ b P.++ c P.++ ")"
  blockStem   _ _ _ w _ _ = w
  blockUnpair _ _ b = "." P.++ b
  compsBR     _ b r = b P.++ (P.replicate (VU.length r) '.')
  compsBC     _ b c = b P.++ c
  structW     _ w   = w
  structCS    _ _ w = "." P.++ w
  structWS    _ w s = w P.++ s
  structOpen  _ r   = P.replicate (VU.length r) '.'
  h = return . id

type CombSignature m e b = Signature m (e, m (SM.Stream m b)) (SM.Stream m b)


(<**)
  :: (Monad m, Eq b, Eq e) -- , Show e, Show (m [b]))
  => Signature m e e
  -> Signature m b (SM.Stream m b)
  -> CombSignature m e b
(<**) f s = (hairpin,interior,multi,blockStem,blockUnpair,compsBR,compsBC,structW,structCS,structWS,structOpen,h) where
  (hairpinF,interiorF,multiF,blockStemF,blockUnpairF,compsBRF,compsBCF,structWF,structCSF,structWSF,structOpenF,hF) = f
  (hairpinS,interiorS,multiS,blockStemS,blockUnpairS,compsBRS,compsBCS,structWs,structCSS,structWSS,structOpenS,hS) = s
  
  xs >>>= f = xs >>= return . SM.map f
  ccm2 xs ys f = xs >>= \xx -> ys >>= \yy -> return $ SM.concatMap (\x -> SM.map (\y -> f x y) yy) xx

  hairpin ener l lp xs rp r = (hairpinF ener l lp xs rp r, return $ SM.singleton $ hairpinS ener l lp xs rp r)
  interior ener l ls li (wF,wS) ri rs r = (interiorF ener l ls li wF ri rs r, wS >>>= \w -> interiorS ener l ls li w ri rs r)
  multi ener l li (bF,bS) (cF,cS) ri r = (multiF ener l li bF cF ri r, ccm2 bS cS $ \b c -> multiS ener l li b c ri r)
  blockStem ener lo l (sF,sS) r ro = (blockStemF ener lo l sF r ro, sS >>>= \s -> blockStemS ener lo l s r ro)
  blockUnpair ener c (bF,bS) = (blockUnpairF ener c bF, bS >>>= \s -> blockUnpairS ener c s)
  compsBR ener (bF,bS) reg = (compsBRF ener bF reg, bS >>>= \s -> compsBRS ener s reg)
  compsBC ener (bF,bS) (cF,cS) = (compsBCF ener bF cF, ccm2 bS cS $ \b c -> compsBCS ener b c)
  structW ener (wF,wS) = (structWF ener wF, wS >>>= \w -> structWs ener w)
  structCS ener c (wF,wS) = (structCSF ener c wF, wS >>>= \w -> structCSS ener c w)
  structWS ener (wF,wS) (sF,sS) = (structWSF ener wF sF, ccm2 wS sS $ \w s -> structWSS ener w s)
  structOpen ener r = (structOpenF ener r, return . SM.singleton $ structOpenS ener r)
  h xs = do
    hfs <- hF $ SM.map P.fst xs
    let phfs = SM.concatMapM P.snd . SM.filter ((hfs==) . P.fst) $ xs
    hS phfs

rnaFoldConstrained :: Vienna2004 -> Primary -> D1Secondary -> (Deka,[String])
rnaFoldConstrained ener inp s = (struct ! (Z:.subword 0 n), bt) where
  (_,Z:.Subword (_:.n)) = bounds weak
  len = VU.length inp
  (weak,block,comps,struct) = unsafePerformIO (rnaFoldFill ener (Just s) inp)
  bt = backtrack ener (Just s) inp (weak,block,comps,struct)
{-# NOINLINE rnaFoldConstrained #-}

rnaFold :: Vienna2004 -> Primary -> (Deka,[String])
rnaFold ener inp = (struct ! (Z:.subword 0 n), bt) where
  (_,Z:.Subword (_:.n)) = bounds weak
  len = VU.length inp
  (weak,block,comps,struct) = unsafePerformIO (rnaFoldFill ener Nothing inp)
  bt = backtrack ener Nothing inp (weak,block,comps,struct)
{-# NOINLINE rnaFold #-}

rnaFoldFill :: Vienna2004 -> Maybe (D1Secondary) -> Primary -> IO (PA.Unboxed (Z:.Subword) Deka, PA.Unboxed (Z:.Subword) Deka, PA.Unboxed (Z:.Subword) Deka, PA.Unboxed (Z:.Subword) Deka)
rnaFoldFill !ener !cs !inp = do
  let n = VU.length inp
  !weak'  <- newWithM (Z:.subword 0 0) (Z:.subword 0 n) huge
  !block' <- newWithM (Z:.subword 0 0) (Z:.subword 0 n) huge
  !comps' <- newWithM (Z:.subword 0 0) (Z:.subword 0 n) huge
  !struc' <- newWithM (Z:.subword 0 0) (Z:.subword 0 n) 0
  fillTables $ gRNAfold ener mfe (mTblSw NonEmptyT weak') (mTblSw NonEmptyT block') (mTblSw NonEmptyT comps') (mTblSw NonEmptyT struc') cs inp
  weakF  <- freeze weak'
  blockF <- freeze block'
  compsF <- freeze comps'
  strucF <- freeze struc'
  return (weakF,blockF,compsF,strucF)
{-# NOINLINE rnaFoldFill #-}

fillTables (MTbl _ weak, weakF, MTbl _ block, blockF, MTbl _ comps, compsF, MTbl _ struc, strucF) = do
  let (_,Z:.Subword (0:.n)) = boundsM weak
  forM_ [n,n-1..0] $ \i -> forM_ [i..n] $ \j -> do
    weakF (subword i j) >>= writeM weak (Z:.subword i j)
    blockF (subword i j) >>= writeM block (Z:.subword i j)
    compsF (subword i j) >>= writeM comps (Z:.subword i j)
    strucF (subword i j) >>= writeM struc (Z:.subword i j)
{-# INLINE fillTables #-}

-- * backtracking

backtrack ener cs (inp :: Primary) (weak :: PA.Unboxed (Z:.Subword) Deka, block :: PA.Unboxed (Z:.Subword) Deka, comps :: PA.Unboxed (Z:.Subword) Deka, struct :: PA.Unboxed (Z:.Subword) Deka) = unId . SM.toList . unId $ sF $ subword 0 n where
  n = VU.length inp
  w :: SwBtTbl Id Deka String
  w = btTbl NonEmptyT weak   (wF :: Subword -> Id (SM.Stream Id String))
  b :: SwBtTbl Id Deka String
  b = btTbl NonEmptyT block  (bF :: Subword -> Id (SM.Stream Id String))
  c :: SwBtTbl Id Deka String
  c = btTbl NonEmptyT comps  (cF :: Subword -> Id (SM.Stream Id String))
  s :: SwBtTbl Id Deka String
  s = btTbl NonEmptyT struct (sF :: Subword -> Id (SM.Stream Id String))
  (_,wF,_,bF,_,cF,_,sF) = gRNAfold ener (mfe <** pretty) w b c s cs inp
{-# INLINE backtrack #-}

