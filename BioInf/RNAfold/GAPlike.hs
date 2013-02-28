{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE EmptyDataDecls #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE PackageImports #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}


module BioInf.RNAfold.GAPlike where

import Data.Vector.Unboxed (Vector(..))
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Fusion.Stream.Monadic as S
import Data.Vector.Fusion.Stream.Size
import Control.Lens
import Data.Array.Repa.Index
import Control.Monad.ST
import Control.Monad.Primitive
import Control.Monad

import ADP.Fusion.GAPlike as ADP
import Data.PrimitiveArray hiding ((!))
import qualified Data.PrimitiveArray as PA
import qualified Data.PrimitiveArray.Zero as PA

import Biobase.Primary
import Biobase.Turner
import Biobase.Vienna

import qualified Biobase.Turner.Import as TI



infixl 8 !
(!) = (PA.!)

-- | The signature for RNAfold. This one is much more complicated than for
-- Nussinov78. We have to deal with many ways on how to score structures and
-- need ways to "look into" the inner part of the sequence.

data Signature m a r = Signature
  { hairpin :: Nuc -> Nuc -> Vector Nuc -> Nuc -> Nuc -> a
  , stacked :: Nuc -> Nuc -> Vector Nuc -> a -> Vector Nuc -> Nuc -> Nuc -> a
  , bulgel  :: Nuc -> Nuc -> Vector Nuc -> a -> Nuc -> Nuc -> a
  , bulger  :: Nuc -> Nuc -> a -> Vector Nuc -> Nuc -> Nuc -> a
  , iloop1N :: Nuc -> Nuc -> a -> Vector Nuc -> Nuc -> Nuc -> a
  , iloopN1 :: Nuc -> Nuc -> Vector Nuc -> a -> Nuc -> Nuc -> a
  , iloop   :: Nuc -> Nuc -> Vector Nuc -> Nuc -> Nuc -> a -> Nuc -> Nuc -> Vector Nuc -> Nuc -> Nuc -> a
  , mloop   :: Nuc -> Nuc -> a -> a -> Nuc -> Nuc -> a
  , blockStem :: Nuc -> Nuc -> a -> Nuc -> Nuc -> a
  , blockBlock :: Nuc -> a -> a
  , compBlock :: a -> Vector Nuc -> a
  , compBlCmp :: a -> a -> a
  , iD :: a -> a
  , bStruct :: Nuc -> a -> a
  , wStruct :: a -> a -> a
  , h :: S.Stream m a -> m r
  }

-- | Minimum-free energy algebra.

mfe :: (Monad m) => Turner2004Model Deka -> Signature m Deka Deka
mfe t = Signature
  { hairpin = \ lo li xs ri ro ->
      let
        len = VU.length xs
        zl  = Z:.len
        p   = Z:.lo:.ro:.li:.ri
      in if
        | len <  3 -> 999999
        | len == 3 -> t ^. termAU + t^.hairpinL ! zl
        | len > 30 -> 999999 -- TODO
        | otherwise -> t^.hairpinL ! zl + t^.hairpinMM ! p
  , stacked = \ !lo !li !xs !z !ys !ri !ro ->
      let
        lxs = VU.length xs
        lys = VU.length ys
        p   = Z:.lo:.ro:.ri:.li
      in if
        | lxs == 0 && lys == 0 -> z + t^.stack ! p
        -- TODO continue with other small loops
        | otherwise -> 999999
  , bulgel = \ lo li xs z ri ro ->
      let
        len = VU.length xs
        l = Z:.len
      in if
        | len <= 30 -> z + t^.bulgeL ! l -- TODO terminal AU inner outer
        | otherwise -> 999999
  , bulger = \ lo li z xs ri ro ->
      let
        len = VU.length xs
        l = Z:.len
      in if
        | len <= 30 -> z + t^.bulgeL ! l
        | otherwise -> 999999
  , iloop1N = \ lo li z xs ri ro ->
      let
        len = VU.length xs
        l = Z:.len
      in if
        | otherwise -> 999999
  , iloopN1 = \ lo li xs z ri ro ->
      let
      in if
        | otherwise -> 999999
  , iloop = \ lo li xs ilo ili z iri iro ys ri ro ->
      let
      in if
        | otherwise -> 999999
  , mloop = \ lo li x y ri ro ->
      let
      in if
        | otherwise -> 999999
  , blockStem = \ lo li z ri ro ->
      let
      in if
        | otherwise -> 999999
  , blockBlock = \ b z ->
      let
      in if
        | otherwise -> 999999
  , compBlock = \ z xs ->
      let
      in if
        | otherwise -> 999999
  , compBlCmp = \ x y ->
      let
      in if
        | otherwise -> 999999
  , iD = \ x -> x
  , bStruct = \ b z ->
      let
      in if
        | otherwise -> 999999
  , wStruct = \ x y ->
      let
      in if
        | otherwise -> 999999
  , h = S.foldl' min 999999
  }
{-# INLINE mfe #-}

gRNAfold Signature{..} (w',bl',c',s') (b, pl, pr, region, small) =
  ( w', ( hairpin <<< b % pr % region % pl % b |||
          stacked <<< b % pr % small % w % small % pl % b |||
          bulgel  <<< b % pr % region % w % pl % b |||
          bulger  <<< b % pr % w % region % pl % b |||
          iloop1N <<< b % b % w % region % pl % b |||
          iloopN1 <<< b % pr % region % w % b % b |||
          iloop   <<< b % pr % region % pl % b % w % b % pr % region % pl % b |||
          mloop   <<< b % pr % bl % c % pl % b ... h
        )
  , bl',  ( blockStem  <<< pl % b % w % b % pr |||
            blockBlock <<< b % bl ... h
          )
  , c',   ( compBlock <<< bl % region |||
            compBlCmp <<< bl % c |||
            iD        <<< bl ... h
          )
  , s',   ( iD      <<< w |||
            bStruct <<< b % s |||
            wStruct <<< w % s ... h
          )
  ) where
      w = mtblN w'
      bl = mtblN bl'
      c = mtblN c'
      s = mtblN s'
{-# INLINE gRNAfold #-}

-- |

test = (f <<< r % r ... S.toList) (0,3) where
  f = (,)
  p = Chr inp
  pl = Peek L inp
  pr = Peek R inp
  r = Regioned 0 5 inp
  inp = VU.fromList [0 .. 99]

data LR = L | R
  deriving (Eq,Show)

data Peek e = Peek !LR !(VU.Vector e)

instance Build (Peek e)

instance (StreamElement x) => StreamElement (x:.Peek e) where
  data StreamElm (x:.Peek e) = SePeek !(StreamElm x) !Int !e
  type StreamTopIdx (x:.Peek e) = Int
  type StreamArg (x:.Peek e) = StreamArg x :. e
  getTopIdx (SePeek _ k _) = k
  getArg (SePeek x _ e) = getArg x :. e
  {-# INLINE getTopIdx #-}
  {-# INLINE getArg #-}

-- TODO I think, we can rewrite both versions to use S.map instead of S.flatten.

instance (Monad m, MkStream m x, StreamElement x, StreamTopIdx x ~ Int, VU.Unbox e) => MkStream m (x:.Peek e) where
  mkStream (x:.Peek lr es) (i,j) = S.flatten mk step Unknown $ mkStream x (i,j) where
    mk :: StreamElm x -> m (StreamElm x, Int)
    mk !x = return (x, getTopIdx x)
    step :: (StreamElm x, Int) -> m (S.Step (StreamElm x, Int) (StreamElm (x:.Peek e)))
    step (!x,!k)
      | i<j && k <= j && lr==L = return $ S.Yield (SePeek x k (VU.unsafeIndex es (j-1))) (x,j+1)
      | i<j && k <= j && lr==R && j+1 < VU.length es = return $ S.Yield (SePeek x k (VU.unsafeIndex es j)) (x,j+1)
      | otherwise = return S.Done
  {-# INLINE mkStream #-}
  mkStreamInner (x:.Peek lr es) (i,j) = S.flatten mk step Unknown $ mkStreamInner x (i,j) where
    mk :: StreamElm x -> m (StreamElm x, Int)
    mk !x = return (x, getTopIdx x)
    step :: (StreamElm x, Int) -> m (S.Step (StreamElm x, Int) (StreamElm (x:.Peek e)))
    step (!x,!k)
      | k < j && k>0 && lr==L = return $ S.Yield (SePeek x k (VU.unsafeIndex es (k-1))) (x,j+1)
      | k < j && lr==R = return $ S.Yield (SePeek x k (VU.unsafeIndex es k)) (x,j+1)
      | otherwise      = return $ S.Done
  {-# INLINE mkStreamInner #-}

-- |

data Region c = Region (VU.Vector c)

instance Build (Region c)

instance (StreamElement x) => StreamElement (x:.Region e) where
  data StreamElm (x:.Region e) = SeRegion !(StreamElm x) !Int !(VU.Vector e)
  type StreamTopIdx (x:.Region e) = Int
  type StreamArg (x:.Region e) = StreamArg x :. (VU.Vector e)
  getTopIdx (SeRegion _ k _) = k
  getArg (SeRegion x _ e) = getArg x :. e
  {-# INLINE getTopIdx #-}
  {-# INLINE getArg #-}

instance (Monad m, MkStream m x, StreamElement x, StreamTopIdx x ~ Int, VU.Unbox e) => MkStream m (x:.Region e) where
  mkStream (x:.Region e) (i,j) = S.map step $ mkStreamInner x (i,j) where
    step :: StreamElm x -> StreamElm (x:.Region e)
    step !x = let k = getTopIdx x in SeRegion x j (VU.unsafeSlice k (j-k) e)
  -- | The inner stream will, in each step, check if the current subword [k,l]
  -- (forall l>=k) is valid and terminate the stream once l>j.
  mkStreamInner (x:.Region e) (i,j) = S.flatten mk step Unknown $ mkStreamInner x (i,j) where
    mk :: StreamElm x -> m (StreamElm x, Int)
    mk !x = return (x, getTopIdx x)
    step :: (StreamElm x, Int) -> m (S.Step (StreamElm x, Int) (StreamElm (x:.Region e)))
    step (!x,!l)
      | l<=j      = return $ S.Yield (SeRegion x l (VU.unsafeSlice k (l-k) e)) (x,l+1)
      | otherwise = return $ S.Done
      where k = getTopIdx x
  {-# INLINE mkStream #-}
  {-# INLINE mkStreamInner #-}

-- |

data Regioned c = Regioned Int Int (VU.Vector c)

instance Build (Regioned c)

instance (StreamElement x) => StreamElement (x:.Regioned e) where
  data StreamElm (x:.Regioned e) = SeRegioned !(StreamElm x) !Int !(VU.Vector e)
  type StreamTopIdx (x:.Regioned e) = Int
  type StreamArg (x:.Regioned e) = StreamArg x :. (VU.Vector e)
  getTopIdx (SeRegioned _ k _) = k
  getArg (SeRegioned x _ e) = getArg x :. e
  {-# INLINE getTopIdx #-}
  {-# INLINE getArg #-}

instance (Monad m, MkStream m x, StreamElement x, StreamTopIdx x ~ Int, VU.Unbox e) => MkStream m (x:.Regioned e) where
    {-
  mkStream (x:.Regioned _ _ e) (i,j) = S.map step $ mkStreamInner x (i,j) where
    step :: StreamElm x -> StreamElm (x:.Regioned e)
    step !x = let k = getTopIdx x in SeRegioned x j (VU.unsafeSlice k (j-k) e)
      -}
  mkStream (x:.Regioned minL maxL e) (i,j) = S.flatten mk step Unknown $ mkStreamInner x (i,j-minL) where
    mk :: StreamElm x -> m (StreamElm x, Int)
    mk !x = return (x, getTopIdx x)
    step :: (StreamElm x, Int) -> m (S.Step (StreamElm x, Int) (StreamElm (x:.Regioned e)))
    step (!x,!k)
      | k<=j && j-k > maxL  = return $ S.Done -- S.Skip (x,k+1)
      | k<=j && j-k >= minL  = return $ S.Yield (SeRegioned x j (VU.unsafeSlice k (j-k) e)) (x,j+1)
      | otherwise   = return $ S.Done
  mkStreamInner (x:.Regioned minL maxL e) (i,j) = S.flatten mk step Unknown $ mkStreamInner x (i,j-minL) where
    mk :: StreamElm x -> m (StreamElm x, Int)
    mk !x = return (x, getTopIdx x)
    step :: (StreamElm x, Int) -> m (S.Step (StreamElm x, Int) (StreamElm (x:.Regioned e)))
    step (!x,!l)
      | l-k <= minL       = return $ S.Done -- S.Skip (x,l+1)
      | l<=j && l-k<=maxL = return $ S.Yield (SeRegioned x l (VU.unsafeSlice k (l-k) e)) (x,l+1)
      | otherwise         = return $ S.Done
      where k = getTopIdx x
  {-# INLINE mkStream #-}
  {-# INLINE mkStreamInner #-}

-- |

run inp = do
  tm <- fmap viennaModel $ TI.fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
  print $ rnafold tm inp

rnafold tm inp = (s ! (Z:.0:.n)) where
  (_,Z:._:.n) = bounds s
  len = VU.length vinp
  vinp = mkPrimary inp
  (w,b,c,s) = runST (rnafoldFill tm vinp)
{-# NOINLINE rnafold #-}

-- |

rnafoldFill :: forall s . Turner2004Model Deka -> Primary -> ST s (PA.U DIM2 Deka, PA.U DIM2 Deka, PA.U DIM2 Deka, PA.U DIM2 Deka)
rnafoldFill tm inp = do
  let n = VU.length inp
  w <- fromAssocsM (Z:.0:.0) (Z:.n:.n) 0 []
  b <- fromAssocsM (Z:.0:.0) (Z:.n:.n) 0 []
  c <- fromAssocsM (Z:.0:.0) (Z:.n:.n) 0 []
  s <- fromAssocsM (Z:.0:.0) (Z:.n:.n) 0 []
  let p = Chr inp
  let pl = Peek L inp
  let pr = Peek R inp
  let region = Region inp
  let r30 = Regioned 0 30 inp
  let small = Region inp -- Regioned 0 4 inp
  fillTables $ gRNAfold (mfe tm) (w,b,c,s) (p,pl,pr,r30,small)
  w' <- freeze w
  b' <- freeze b
  c' <- freeze c
  s' <- freeze s
  return (w',b',c',s')
{-# NOINLINE rnafoldFill #-}

-- |

fillTables :: PrimMonad m =>
  ( PA.MU m DIM2 Deka, (Int,Int) -> m Deka
  , PA.MU m DIM2 Deka, (Int,Int) -> m Deka
  , PA.MU m DIM2 Deka, (Int,Int) -> m Deka
  , PA.MU m DIM2 Deka, (Int,Int) -> m Deka
  ) -> m ()
fillTables
  ( w, wf
  , b, bf
  , c, cf
  , s, sf
  ) = do
  let (_,Z:.n:._) = boundsM w
  forM_ [n,n-1..0] $ \i -> forM_ [i..n] $ \j -> do
    wv <- wf (i,j)
    wv `seq` writeM w (Z:.i:.j) wv
    bv <- bf (i,j)
    bv `seq` writeM b (Z:.i:.j) bv
    cv <- cf (i,j)
    cv `seq` writeM c (Z:.i:.j) cv
    sv <- sf (i,j)
    sv `seq` writeM s (Z:.i:.j) sv
{-# INLINE fillTables #-}

