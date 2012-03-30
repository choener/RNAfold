{-# LANGUAGE PackageImports #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE RecordWildCards #-}

-- | A set of energy functions that are modelled after the ViennaRNA package,
-- Version 2 (with d=2).
--
-- As part of the design, we could have (i) either continued giving many
-- parameters or (ii) have fewer parameters that require input of the
-- (Primary,Index,Index) type. Since the compilation speed of the grammars
-- using these functions depends on the number of arguments (in case (i)
-- compilation takes minutes!), this approach has benefits for testing the
-- fusion library.
--
-- This means compilation is a lot faster, but runtime is not @2.8x@ slower but
-- @3.5x@ slower.

module BioInf.RNAfold.Energy where

import Control.Exception (assert)
import Data.Strict.Tuple hiding (fst,snd)
import qualified Data.Vector.Fusion.Stream.Monadic as S

import Biobase.Primary
import Biobase.Secondary.Vienna
import Biobase.Vienna
import Data.PrimitiveArray
import "PrimitiveArray" Data.Array.Repa.Index

import Debug.Trace



-- | Hairpin structures. Hairpins with less than 3 unpaired nucleotides are
-- forbidden.
--
-- NOTE @(xs,i,j)@ is indeed *only* the unpaired stretch. Hence, the length is
-- @j-i+1@, as given.
--
-- TODO Activate tabulated hairpin structures.

hairpinF :: Vienna2004 -> (Nuc :!: Nuc) -> (Primary :!: Int :!: Int) -> (Nuc :!: Nuc) -> Int
hairpinF Vienna2004{..} (l :!: lp) (xs :!: i :!: j) (rp :!: r) = hairpinType where
  {-# INLINE hairpinType #-}
  hairpinType
    -- TODO use sliceEq for tabulated hairpins
    {-
    | len <= 6, Just v <- find tabulated hairpin
    -}
    | len <  3  = 999999
    | len == 3  = hairpinL!(Z:.len) + tAU
    | len > 31  = hairpinL!(Z:.30)  + hairpinMM!(Z:.p:.lp:.rp) + llp
    | otherwise = hairpinL!(Z:.len) + hairpinMM!(Z:.p:.lp:.rp)
  p   = mkViennaPair (l,r)
  len = j-i+1
  tAU = if p/=vpCG && p/=vpGC then termAU else 0
  llp = floor $ 108.856 * log (fromIntegral len / 30)
{-# INLINE hairpinF #-}

-- | Tiny loops are small interior loops. This includes canonical stacks
-- without any unpaired nucleotides, and small, tabulated interior loops.

tinyloopF :: Vienna2004 -> (Primary :!: Int :!: Int) -> Int -> (Primary :!: Int :!: Int) -> Int
tinyloopF Vienna2004{..} (ls :!: li :!: lj) w (rs :!: ri :!: rj) = assert (assertL && assertR) $ loopType where
  assertL = let (_,Z:.n) = bounds ls in li>=0 && lj<=n
  assertR = let (_,Z:.n) = bounds rs in ri>=0 && rj<=n
  {-# INLINE loopType #-}
  loopType
    | dl==0 || dr==0 = error $ "bug in tinyloop: " ++ show (li,lj,[ls!(Z:.k) | k<-[li..lj]],ri,rj,[rs!(Z:.l) | l<-[ri..rj]])
    -- normal stack
    | dl==1 && dr==1 = w+stack!(Z:.op:.ip)
    -- one intervening unpaired nucleotide doesn't break the stack
    | dl==1 && dr==2 = w+stack!(Z:.op:.ip)+bulgeL!(Z:.1)
    | dl==2 && dr==1 = w+stack!(Z:.op:.ip)+bulgeL!(Z:.1)
    -- 1x1 symmetric interior loop
    | dl==2 && dr==2 = w+iloop1x1!(Z:.op:.ip:.bI:.bJ)
    -- 1x2 interior loop
    | dl==2 && dr==3 = w+iloop2x1!(Z:.op:.ip:.bI:.bL:.bJ)
    -- 2x1 interior loop
    | dl==3 && dr==2 = w+iloop2x1!(Z:.op:.ip:.bJ:.bI:.bK)
    -- 2x2 interior loop
    | dl==3 && dr==3 = w+iloop2x2!(Z:.op:.ip:.bI:.bK:.bL:.bJ)
    -- 2x3 interior loops with specialized mismatches
    |  dl==3 && dr==4
    || dl==4 && dr==3 = w+iloop2x3MM!(Z:.op:.bI:.bJ) + iloop2x3MM!(Z:.ip:.bL:.bK) + iloopL!(Z:.5) + ninio
    | otherwise      = 999999 -- overlaps with big interior loops
  op = mkViennaPair (ls!(Z:.li),rs!(Z:.rj))
  ip = mkViennaPair (rs!(Z:.ri),ls!(Z:.lj))
  bI = ls!(Z:.li+1)
  bJ = rs!(Z:.rj-1)
  bK = ls!(Z:.lj-1)
  bL = rs!(Z:.ri+1)
  dl = lj-li
  dr = rj-ri
{-# INLINE tinyloopF #-}

-- | A left bulge @(....[[...]])@ which four unpaired nucleotides in the bulge.
-- the left bulge @ls@ will be given six nucleotides (note, @ls@ is the
-- complete input, use @li@ and @lj@ as the first and last included nucleotide
-- index), the two outer ones being for the outer and inner loop. On the right,
-- we have @rp@ and @r@ which are nucleotides. @ls!(Z:.li)@ and @r@ form the
-- outer Vienna pair. @rp@ and @ls!(Z:.lj)@ form the inner pair.

bulgeLF :: Vienna2004 -> (Primary :!: Int :!: Int) -> Int -> (Nuc :!: Nuc) -> Int
bulgeLF Vienna2004{..} (ls :!: li :!: lj) w (rp :!: r) = assert (lj-li<=30) $ w + tAUlr + lenE + tAUrpcp where
  tAUlr = terminalAU termAU lr
  tAUrpcp = terminalAU termAU rpcp
  lr = mkViennaPair (ls!(Z:.li),r)
  rpcp = mkViennaPair (rp,ls!(Z:.lj))
  lenE = bulgeL!(Z:.lj-li-1)
{-# INLINE bulgeLF #-}

-- | A right bulge @([[...]]....)@. See 'bulgeLF' for how this works.

bulgeRF :: Vienna2004 -> (Nuc :!: Nuc) -> Int -> (Primary :!: Int :!: Int) -> Int
bulgeRF Vienna2004{..} (l :!: lp) w (rs :!: ri :!: rj) = assert (rj-ri<=30) $ w + tAUlr + lenE + tAUcplp where
  tAUlr = terminalAU termAU lr
  tAUcplp = terminalAU termAU cplp
  lr = mkViennaPair (l,rs!(Z:.rj))
  cplp = mkViennaPair (rs!(Z:.ri),lp)
  lenE = bulgeL!(Z:.rj-ri-1)
{-# iNLINE bulgeRF #-}

-- | An interior loop with @N@ unpaired nucleotides to the left and @1@
-- unpaired nucleotide to the right. The regions @ls@ and @rs@ each have 2
-- nucleotides more than are unpaired. These first and last nucleotides form
-- the last paired or first pairs in the stacks around the loop.

iloopN1F :: Vienna2004 -> (Primary :!: Int :!: Int) -> Int -> (Primary :!: Int :!: Int) -> Int
iloopN1F Vienna2004{..} (ls:!:li:!:lj) w (rs:!:ri:!:rj)
  = assert (lj-li-1 <=29 && rj-ri == 2)
  $ w + outerMM + lenE + innerMM + owninio where
    poLR = mkViennaPair (ls!(Z:.li), rs!(Z:.rj))
    piRL = mkViennaPair (rs!(Z:.ri), ls!(Z:.lj))
    outerMM = iloop1xnMM!(Z:.poLR:.ls!(Z:.li+1):.rs!(Z:.ri+1))
    innerMM = iloop1xnMM!(Z:.piRL:.rs!(Z:.ri+1):.ls!(Z:.lj-1))
    lenE = iloopL!(Z:.(lL+1))
    owninio = min maxNinio (ninio * (lL-1))
    lL = lj-li-1
    lR = rj-ri-1
{-# INLINE iloopN1F #-}

-- | 1xN interior loops.

iloop1NF :: Vienna2004 -> (Primary :!: Int :!: Int) -> Int -> (Primary :!: Int :!: Int) -> Int
iloop1NF Vienna2004{..} (ls :!: li :!: lj) w (rs :!: ri :!: rj)
  = assert (lj-li == 2 && rj-ri-1 <=29)
  $ w + outerMM + lenE + innerMM + owninio where
    poLR = mkViennaPair (ls!(Z:.li), rs!(Z:.rj))
    piRL = mkViennaPair (rs!(Z:.ri), ls!(Z:.lj))
    outerMM = iloop1xnMM!(Z:.poLR:.ls!(Z:.li+1):.rs!(Z:.rj-1))
    innerMM = iloop1xnMM!(Z:.piRL:.rs!(Z:.ri+1):.ls!(Z:.li+1))
    lenE = iloopL!(Z:.(lR+1))
    owninio = min maxNinio (ninio * (lR-1))
    lL = lj-li-1
    lR = rj-ri-1
{-# INLINE iloop1NF #-}

iloopIF :: Vienna2004 -> (Primary :!: Int :!: Int) -> Int -> (Primary :!: Int :!: Int) -> Int
iloopIF Vienna2004{..} (ls :!: li :!: lj) w (rs :!: ri :!: rj) = w + iloopMM!(Z:.p:.bR:.bL) + iloopL!(Z:.len) + owninio where
  p = mkViennaPair (rs!(Z:.ri), ls!(Z:.lj))
  bL = ls!(Z:.lj-1)
  bR = rs!(Z:.ri+1)
  len = (lj-li)+(rj-ri)
  owninio = min maxNinio (ninio * (abs $ (lj-li) - (rj-ri)))
{-# INLINE iloopIF #-}

-- |

iloopOF :: Vienna2004 -> (Nuc :!: Nuc) -> Int -> (Nuc :!: Nuc) -> Int
iloopOF Vienna2004{..} (l :!: lp) iloopif (rp :!: r) = {- CORE "iloopOF" -} iloopif + iloopMM!(Z:.op:.lp:.rp)
  where op = mkViennaPair (l,r)
{-# INLINE iloopOF #-}

-- |

multiIF :: Vienna2004 -> Int -> Int -> Int
multiIF Vienna2004{..} b c = b+c
{-# INLINE multiIF #-}

-- |

multiOF :: Vienna2004 -> (Nuc :!: Nuc) -> Int -> (Nuc :!: Nuc) -> Int
multiOF Vienna2004{..} (l :!: lp) multiif (rp :!: r) = multiif + multiMM!(Z:.p:.lp:.rp) + multiHelix + multiOffset + terminalAU termAU p
  where p = mkViennaPair (l,r)
{-# INLINE multiOF #-}

-- |

regionStemF :: Vienna2004 -> Nuc -> Int -> Int
regionStemF Vienna2004{..} _ w = w + multiNuc
{-# INLINE regionStemF #-}

-- |

justStemF :: Vienna2004 -> (Nuc :!: Nuc) -> Int -> (Nuc :!: Nuc) -> Int
justStemF Vienna2004{..} (l :!: lp) w (rp :!: r) = w + multiMM!(Z:.p:.l:.r)
  where p = mkViennaPair (lp,rp)
{-# INLINE justStemF #-}

-- |

bsF :: Int -> Int -> Int
bsF b reg = b
{-# INLINE bsF #-}

-- |

rSF :: Nuc -> Int -> Int
rSF nuc w = w
{-# INLINE rSF #-}

-- |

bcF :: Int -> Int -> Int
bcF b c = b + c
{-# INLINE bcF #-}

-- |

ssF :: Int -> Int
ssF reg = 0
{-# INLINE ssF #-}

-- |

cmF :: Int -> Int -> Int
cmF w s = w+s
{-# INLINE cmF #-}

-- |

nilF :: Bool -> Int
nilF e = if e then 0 else 999999
{-# INLINE nilF #-}

-- |

iD = id
{-# INLINE iD #-}

-- |

h :: Monad m => S.Stream m Int -> m Int
h = S.foldl' min (999999::Int)
{-# INLINE h #-}

-- |

terminalAU termAU p = if p/=vpCG && p/=vpGC then termAU else 0
{-# INLINE terminalAU #-}

