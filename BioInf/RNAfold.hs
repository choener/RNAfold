{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RecordWildCards #-}

module BioInf.RNAfold where

import Control.Monad
import Control.Monad.Primitive
import Control.Monad.ST
import Data.Array.Repa.Index
import qualified Data.Vector.Fusion.Stream.Monadic as S
import qualified Data.Vector.Fusion.Stream as P
import qualified Data.Vector.Unboxed as VU
import Control.Monad.State.Lazy
import Control.Arrow (first,second,(***))
import qualified Data.List as L
import Control.Exception (assert)
import Data.Strict.Tuple hiding (fst,snd)

import Biobase.Primary
import Biobase.Secondary.Vienna

import Data.PrimitiveArray
import Data.PrimitiveArray.Unboxed.Zero

import ADP.Fusion.Monadic
import ADP.Fusion.Monadic.Internal

import Debug.Trace
import Text.Printf
import GHC.Exts

import Biobase.Primary
import Biobase.Secondary.Vienna
import Biobase.Vienna
import Biobase.Vienna.Default

import BioInf.RNAfold.Combinators
import BioInf.RNAfold.Energy
import BioInf.RNAfold.Library



-- |

testRNAfold :: String -> (Int,[String])
testRNAfold inp' = struct `seq` bt `seq` (struct!(Z:.0:.n),bt) where
  (_,Z:._:.n) = bounds struct
  ener = fst turnerRNA2004
  inp = mkPrimary inp'
  tbls@(weak,block,comps,struct) = runST (rnafold ener inp)
  bt = btRNAfold ener inp tbls
{-# NOINLINE testRNAfold #-}

-- |
testInput = "cccacccaaagggaaaaggg"
test = testRNAfold testInput

rnafold :: Vienna2004 -> Primary -> ST s
  ( Arr0 DIM2 Int
  , Arr0 DIM2 Int
  , Arr0 DIM2 Int
  , Arr0 DIM2 Int
  )
rnafold ener inp = do

  let !n = let (_,Z:.l) = bounds inp in l+1
  let base   = base'   inp
      {-# INLINE base #-}
  let baseLr = baseLr' inp
      {-# INLINE baseLr #-}
  let baselR = baselR' inp
      {-# INLINE baselR #-}
  let basepairing = basepairing' inp
      {-# INLINE basepairing #-}
  let stackpairing = stackpairing' inp
      {-# INLINE stackpairing #-}
  let reglen = reglen' inp
      {-# INLINE reglen #-}
  let primary = primary' inp
      {-# INLINE primary #-}
  let primaryPR = primaryPR' inp
      {-# INLINE primaryPR #-}
  let primaryPL = primaryPL' inp
      {-# INLINE primaryPL #-}
  -- this is a bit unfortunate, but otherwise we get type inference problems
  let hS :: S.Stream (ST s) Int -> ScalarM (ST s Int)
      hS = ScalarM . S.foldl' min (999999::Int)
      {-# INLINE hS #-}

  weak   :: MArr0 s DIM2 Int <- fromAssocsM (Z:.0:.0) (Z:.n:.n) 999999 []
  block  :: MArr0 s DIM2 Int <- fromAssocsM (Z:.0:.0) (Z:.n:.n) 999999 []
  comps  :: MArr0 s DIM2 Int <- fromAssocsM (Z:.0:.0) (Z:.n:.n) 999999 []
  struct :: MArr0 s DIM2 Int <- fromAssocsM (Z:.0:.0) (Z:.n:.n) 0 []

  let iif = iloopIF ener <<< primary #~~ weak ~~# primary ... hS
      {-# INLINE [0] iif #-}

  let mif = multiIF ener <<< block +~+ comps ... hS
      {-# INLINE [0] mif #-}

  fillTables
    weak (
      -- multiOF   ener <<< baseLr -~+ (multiIF ener <<< block +~+ comps              ... hS) +~- baselR |||
      multiOF   ener <<< baseLr -~+ mif +~- baselR |||
      -- iloopOF   ener <<< baseLr -~+ (iloopIF ener <<< primary #~~ weak ~~# primary ... hS) +~- baselR |||
      iloopOF   ener <<< baseLr -~+ iif +~- baselR |||
      iloop1NF  ener <<< primary ---~+ weak +~@   primary   |||
      iloopN1F  ener <<< primary   @~+ weak +~--- primary   |||
      bulgeRF   ener <<< baseLr    -~+ weak +~*   primary   |||
      bulgeLF   ener <<< primary   *~+ weak +~-   baselR    |||
      tinyloopF ener <<< primaryPR &~+ weak +~&   primaryPL |||
      hairpinF  ener <<< baseLr -~+ primary +~- baselR ... h `with` basepairing )
    block (
      adjustStream n (justStemF ener <<< baseLr -~+ weak +~- baselR) |||
      regionStemF ener <<< base -~+ block             ... h )
    comps (
      bsF <<< block +~+ reglen |||
      bcF <<< block +~+ comps  |||
      iD  <<< block            ... h )
    struct (
      iD   <<< weak            |||
      rSF  <<< base -~~ struct |||
      cmF  <<< weak +~+ struct |||
      nilF <<< empty           ... h `with` constrained (\(Z:.i:.j) -> j==n) )

  weak'   <- freeze weak
  block'  <- freeze block
  comps'  <- freeze comps
  struct' <- freeze struct

  return
    ( weak'
    , block'
    , comps'
    , struct'
    )
{-# INLINE rnafold #-}

infixl 9 *~+, +~*, &~+, +~&, ---~+, +~@, +~---, @~+

(*~+) = makeLeft_MinRight (3,31) 1
{-# INLINE (*~+) #-}

(+~*) = makeMinLeft_Right 1 (3,31)
{-# INLINE (+~*) #-}

(&~+) = makeLeft_MinRight (1,4) 1
{-# INLINE (&~+) #-}

(+~&) = makeMinLeft_Right 1 (1,4)
{-# INLINE (+~&) #-}

(---~+) = makeLeft_MinRight (3,3) 1
{-# INLINE (---~+) #-}

(+~@) = makeMinLeft_Right 1 (5,31)
{-# INLINE (+~@) #-}

(+~---) = makeMinLeft_Right 1 (3,3)
{-# INLINE (+~---) #-}

(@~+) = makeLeft_MinRight (5,31) 1
{-# INLINE (@~+) #-}

-- * backtracking
--
-- For now, we replicate the grammar as the optimizer is rather fragile

btRNAfold
  :: Vienna2004
  -> Primary
  -> (Arr0 DIM2 Int, Arr0 DIM2 Int, Arr0 DIM2 Int, Arr0 DIM2 Int)
  -> [String]
btRNAfold ener inp (weak,block,comps,struct) = structG (Z:.0:.n) where

  !n = let (_,Z:.l) = bounds inp in l+1
  base   = base'   inp
  baseLr = baseLr' inp
  baselR = baselR' inp
  reglen = reglen' inp
  basepairing = basepairing' inp
  primary = primary' inp
  primaryPR = primaryPR' inp
  primaryPL = primaryPL' inp

  --

  weak' :: DIM2 -> Scalar (Int, [String])
  weak' ij = Scalar (weak!ij, weakG ij)

  block' :: DIM2 -> Scalar (Int, [String])
  block' ij = Scalar (block!ij, blockG ij)

  comps' :: DIM2 -> Scalar (Int, [String])
  comps' ij = Scalar (comps!ij, compsG ij)

  struct' :: DIM2 -> Scalar (Int, [String])
  struct' ij = Scalar (struct!ij, structG ij)

  --

  weakG :: DIM2 -> [String]
  weakG = (
            multiBT ener    <<< baseLr -~+ block'  +~+ comps' +~- baselR             |||
            iloopBT ener    <<< baseLr -~+ primary #~~ weak'  ~~# primary +~- baselR |||
            iloop1NBT ener  <<< primary ---~+ weak'   +~@   primary |||
            iloopN1BT ener  <<< primary @~+   weak'   +~--- primary |||
            bulgeRBT ener   <<< baseLr  -~+   weak'   +~*   primary |||
            bulgeLBT ener   <<< primary *~+   weak'   +~-   baselR  |||
            tinyloopBT ener <<< primaryPR &~+   weak'   +~&   primaryPL |||
            hairpinBT  ener <<< baseLr  -~+   primary +~-   baselR  ..@ (hBT weak) `withBT` basepairing
          )

  blockG :: DIM2 -> [String]
  blockG = (
             adjustStreamBT n (justStemBT   ener <<< baseLr -~+ weak'  +~- baselR) |||
             regionStemBT ener <<< base   -~+  block'             ..@ (hBT block)
           )

  compsG :: DIM2 -> [String]
  compsG = (
             bsBT <<< block' +~+ reglen |||
             bcBT <<< block' +~+ comps' |||
             iDBT <<< block'            ..@ (hBT comps)
           )

  structG :: DIM2 -> [String]
  structG = (
              iDBT  <<< weak'             |||
              rSBT  <<< base  -~~ struct' |||
              cmBT  <<< weak' +~+ struct' |||
              nilBT <<< empty             ..@ (hBT struct `withBT` constrained (\(Z:.i:.j) -> j==n) )
            )

  multiBT ener l (b,bS) (c,cS) r =
    let e = multiOF ener l (multiIF ener b c) r
    in (e, ["("++x++y++")" | x<-bS, y<-cS])
  iloopBT ener lo ls@(_:!:li:!:lj) (w,wS) rs@(_:!:ri:!:rj) ro =
    let e = iloopOF ener lo (iloopIF ener ls w rs) ro
    in (e, L.map (\s -> "("++replicate (lj-li) '.'++"("++s++")"++replicate (rj-ri) '.'++")") wS)
  iloop1NBT ener ls (w,wS) rs@(_:!:ri:!:rj) =
    let e = iloop1NF ener ls w rs
    in (e, L.map (\s -> "(.("++s++")"++replicate (rj-ri-1) '.'++")") wS)
  iloopN1BT ener ls@(_:!:li:!:lj) (w,wS) rs =
    let e = iloopN1F ener ls w rs
    in (e, L.map (\s -> "("++replicate (lj-li-1) '.'++"("++s++").)") wS)
  bulgeRBT ener ls (w,wS) rs@(_:!:ri:!:rj) =
    let e = bulgeRF ener ls w rs
    in (e, L.map (\s -> "("++s++"."++replicate (rj-ri-1) '.'++")") wS)
  bulgeLBT ener ls@(_:!:li:!:lj) (w,wS) rs =
    let e = bulgeLF ener ls w rs
    in (e, L.map (\s -> "("++replicate (lj-li-1) '.'++"."++s++")") wS)
  hairpinBT ener llp reg@(xs:!:i:!:j) rpr =
    let e = hairpinF ener llp reg rpr
    in (e, ["(" ++ replicate (j-i+1) '.' ++ ")"])
  tinyloopBT ener ls@(_:!:li:!:lj) (w,sW) rs@(_:!:ri:!:rj) =
    let e = tinyloopF ener ls w rs
    in (e, L.map (\s -> "("++replicate (lj-li-1) '.'++s++replicate (rj-ri-1) '.'++")") sW)

  regionStemBT ener nc (w,sW) =
    let e = regionStemF ener nc w
    in (e, L.map (\s -> "." ++ s) sW)
  justStemBT ener llp (w,sW) rpr =
    let e = justStemF ener llp w rpr
    in (e, sW)

  bcBT (b,bW) (c,cW) =
    let e = b+c
    in (e,[ x++y | x<-bW, y<-cW ])
  bsBT (b,bW) reg = (b, L.map (++ replicate reg '.') bW)
  iDBT = id

  ssBT len = (ssF len, [replicate len '.'])
  rSBT n (w,wS) =
    let e = rSF n w
    in (e, map ("."++) wS)
  cmBT (w,wS) (s,sS) =
    let e = cmF w s
    in (e, [x++y | x<-wS, y<-sS])
  nilBT b = if b then (nilF b, [""]) else (nilF b, [])

  hBT tbl ij = L.concatMap snd . L.filter ((tbl!ij==).fst) . P.toList


-- * different energy functions (very simplified signature)



-- * Functions that will be part of a bioinformatics DP library



-- **

adjustStream :: Int -> (DIM2 -> S.Stream (ST s) Int) -> DIM2 -> S.Stream (ST s) Int
adjustStream !n sgen (Z:.i:.j)
  | i>0 && j<n = sgen (Z:.i-1:.j+1)
  | otherwise  = S.empty
{-# INLINE adjustStream #-}

adjustStreamBT :: Int -> (DIM2 -> P.Stream elm) -> DIM2 -> P.Stream elm
adjustStreamBT !n sgen (Z:.i:.j)
  | i>0 && j<n = sgen (Z:.i-1:.j+1)
  | otherwise  = P.empty
{-# INLINE adjustStreamBT #-}

infixl 6 .!.
(.!.) stream (h,n) (Z:.i:.j)
  | i>0 && j<n = h $ stream (Z:.i-1:.j+1)
  | otherwise  = h $ S.empty
{-# INLINE (.!.) #-}

-- |

infixl 5 `with`
with xs cond ij = if cond ij then xs ij else return 999999
{-# INLINE with #-}

infixl 5 `withBT`
withBT xs cond ij = if cond ij then xs ij else return []

-- |

fillTables
  :: PrimMonad m
  => MArr0 (PrimState m) DIM2 Int -> (DIM2 -> m Int)
  -> MArr0 (PrimState m) DIM2 Int -> (DIM2 -> m Int)
  -> MArr0 (PrimState m) DIM2 Int -> (DIM2 -> m Int)
  -> MArr0 (PrimState m) DIM2 Int -> (DIM2 -> m Int)
  -> m ()
fillTables aT aF bT bF cT cF dT dF = do
  let (_,Z:.n:._) = boundsM aT
  forM_ [n,n-1 .. 0] $ \i -> forM_ [i..n] $ \j -> do
    let ij = Z:.i:.j
    aF ij >>= writeM aT ij
    bF ij >>= writeM bT ij
    cF ij >>= writeM cT ij
    dF ij >>= writeM dT ij
{-# INLINE fillTables #-}

