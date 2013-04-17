{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RecordWildCards #-}

module BioInf.RNAfold where

import Data.Vector.Fusion.Stream.Monadic as SM
import Data.Array.Repa.Index
import Data.Strict.Tuple
import qualified Data.Vector.Unboxed as VU
import Control.Monad.ST
import Control.Monad
import Prelude as P
import System.IO.Unsafe
import Control.Lens

import Data.Array.Repa.Index.Subword
import Biobase.Primary
import ADP.Fusion
import Biobase.Vienna
import Biobase.Turner
import Data.PrimitiveArray as PA hiding ((!))
import Data.PrimitiveArray.Zero as PA
import qualified Data.PrimitiveArray as PA

infixl 8 !
(!) = (PA.!)
{-# INLINE (!) #-}


type Signature m a r =
  ( Vienna2004 -> Nuc -> Nuc -> Primary -> Nuc -> Nuc -> a
  , Vienna2004 -> Nuc -> Primary -> Nuc -> a -> Nuc -> Primary -> Nuc -> a
  , Vienna2004 -> Nuc -> Nuc -> a -> a -> Nuc -> Nuc -> a
  , Stream m a -> m r
  )

gRNAfold ener (hairpin,interior,multi,h) weak block comps inp =
  ( weak ,
    hairpin  ener <<< c % pr % sr % pl % c            |||
    interior ener <<< c % r % pr % weak % pl % r % c  |||
    multi    ener <<< pl % c % block % comps % c % pr ... h
  ) where c = chr inp
          r = region inp
          pr = peekR inp
          pl = peekL  inp
          sr = sregion 3 30 inp
          {-# INLINE c #-}
          {-# INLINE r #-}
          {-# INLINE pr #-}
          {-# INLINE pl #-}
          {-# INLINE sr #-}
{-# INLINE gRNAfold #-}

mfe :: Monad m => Signature m Deka Deka
mfe = (hairpin,interior,multi,h) where
  hairpin ener l lp xs rp r
      | len <   3 = Deka 888888
      | len ==  3 = (ener^.hairpinL) VU.! len + tAU + Deka 69696969
      | len < 31  = (ener^.hairpinL) VU.! len + ener^.hairpinMM!(Z:.l:.r:.lp:.rp) + Deka 498349834983
      | otherwise = Deka 777777
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
      = w + _iloop2x2 ener ! (Z:.l:.r:.r:.li:.lH:.lL:.rH:.rL)
      | min lls lrs == 2 && max lls lrs == 3
      = w + _iloop2x3MM ener ! (Z:.l:.r:.lH:.lL) + _iloop2x3MM ener ! (Z:.ri:.li:.rL:.rH) + _iloopL ener VU.! 5 + _ninio ener
      | lls==0 && lrs > 1 && lrs <= 30
      = w + tAU + _bulgeL ener VU.! lrs + tUA
      | lrs==0 && lls > 1 && lls <= 30
      = w + tAU + _bulgeL ener VU.! lls + tUA
      | lls==1 && lrs > 2 && lrs <= 30
      = w + error "TODO: 1xn loops"
      | lrs==1 && lls > 2 && lls <= 30
      = w + error "TODO: 1xn loops"
      | lls+lrs <= 30 -- TODO missing support for length constraints ?
      = w + error "TODO: general interior loops"
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
  multi ener lo l b c r ro
    = huge
  h = foldl' min huge
  {-# INLINE hairpin #-}
  {-# INLINE interior #-}
  {-# INLINE multi #-}
  {-# INLINE h #-}
{-# INLINE mfe #-}

huge = Deka 999999
{-# INLINE huge #-}

rnaFold ener inp = (weak ! (Z:.subword 0 n), bt) where
  (_,Z:.Subword (_:.n)) = bounds weak
  len = P.length inp
  vinp = mkPrimary inp
  (weak,block,comps) = unsafePerformIO (rnaFoldFill ener vinp)
  bt = [] :: [String]
{-# NOINLINE rnaFold #-}

rnaFoldFill :: Vienna2004 -> Primary -> IO (PA.Unboxed (Z:.Subword) Deka, PA.Unboxed (Z:.Subword) Deka, PA.Unboxed (Z:.Subword) Deka)
rnaFoldFill !ener !inp = do
  let n = VU.length inp
  !weak'  <- newWithM (Z:.subword 0 0) (Z:.subword 0 n) huge
  !block' <- newWithM (Z:.subword 0 0) (Z:.subword 0 n) huge
  !comps' <- newWithM (Z:.subword 0 0) (Z:.subword 0 n) huge
  fillTables $ gRNAfold ener mfe (MTbl NoEmptyT weak') (MTbl NoEmptyT block') (MTbl NoEmptyT comps') inp
  weakF  <- freeze weak'
  blockF <- freeze block'
  compsF <- freeze comps'
  return (weakF,blockF,compsF)
{-# NOINLINE rnaFoldFill #-}

fillTables ((MTbl _ weak, weakF)) = do
  let (_,Z:.Subword (0:.n)) = boundsM weak
  forM_ [n,n-1..0] $ \i -> forM_ [i..n] $ \j -> do
    !v <- (weakF $ subword i j)
    writeM weak (Z:.subword i j) v
{-# INLINE fillTables #-}

{-

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


-}

