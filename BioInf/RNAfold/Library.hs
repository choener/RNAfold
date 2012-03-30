{-# LANGUAGE TypeOperators #-}

-- | A library of helpers for ADPfusion algorithms.

module BioInf.RNAfold.Library where

import Control.Exception (assert)
import Data.Array.Repa.Index
import Debug.Trace
import qualified Data.Vector.Unboxed as VU
import Data.Strict.Tuple hiding (fst,snd)

import Biobase.Primary
import Biobase.Secondary.Vienna
import Data.PrimitiveArray

import ADP.Fusion.Monadic
import ADP.Fusion.Monadic.Internal



-- |

base' :: Primary -> DIM2 -> (Scalar Nuc)
base' inp (Z:.i:.j) = Scalar $ index inp (Z:.i)
{-# INLINE base' #-}

-- | Nucleotide, second one to the right. The assertion allows a size of one or
-- two to capture special cases of looking outside of the bounds (used by
-- @justStemF@ in RNAfold).

baseLr' :: Primary -> DIM2 -> (Scalar (Nuc :!: Nuc))
baseLr' inp (Z:.i:.j)
  = assert (let (_,Z:.n) = bounds inp in (i+1==j || i+2==j) && i>=0 && i+1<=n)
  . Scalar $ (index inp (Z:.i) :!: index inp (Z:.i+1))
{-# INLINE baseLr' #-}

-- |

baselR' :: Primary -> DIM2 -> (Scalar (Nuc :!: Nuc))
baselR' inp (Z:.i:.j)
  = assert (let (_,Z:.n) = bounds inp in i+1==j && i>0 && i<=n)
  . Scalar $ (index inp (Z:.i-1) :!: index inp (Z:.i))
{-# INLINE baselR' #-}

-- |

region' :: VU.Vector Nuc -> DIM2 -> (Scalar (VU.Vector Nuc))
region' inp (Z:.i:.j)
  = assert (let n = VU.length inp -1 in i<=j && i>=0 && j<=n+1)
  . Scalar $ VU.unsafeSlice i (j-i) inp
{-# INLINE region' #-}

-- | A 'Primary' together with the lowest included nucleotide and the highest
-- included nucleotide.

primary' :: Primary -> DIM2 -> (Scalar (Primary :!: Int :!: Int))
primary' inp (Z:.i:.j)
  = assert (let (_,Z:.u) = bounds inp in i>=0 && j-1<=u && i<=j)
  $ Scalar (inp :!: i :!: j-1)
{-# INLINE primary' #-}

-- | A 'Primary' together with the lowest included nucleotide and the highest
-- included nucleotide.

primaryPR' :: Primary -> DIM2 -> (Scalar (Primary :!: Int :!: Int))
primaryPR' inp (Z:.i:.j)
  = assert (let (_,Z:.u) = bounds inp in i>=0 && j-2<=u && i<=j)
  $ Scalar (inp :!: i :!: j)
{-# INLINE primaryPR' #-}

-- | A 'Primary' together with the lowest included nucleotide and the highest
-- included nucleotide.

primaryPL' :: Primary -> DIM2 -> (Scalar (Primary :!: Int :!: Int))
primaryPL' inp (Z:.i:.j)
  = assert (let (_,Z:.u) = bounds inp in i>0 && j-1<=u && i<=j)
  $ Scalar (inp :!: i-1 :!: j-1)
{-# INLINE primaryPL' #-}

-- | Vector of nucleotides peeking one nucleotide to the left.

regionpl' :: VU.Vector Nuc -> DIM2 -> (Scalar (VU.Vector Nuc))
regionpl' inp (Z:.i:.j)
  = assert (let n = VU.length inp -1 in i-1<=j && i>0 && j<=n+1)
  . Scalar $ VU.unsafeSlice (i-1) (j-i+1) inp
{-# INLINE regionpl' #-}

-- | Vector of nucleotides peeaking one nucleotide to the right.

regionpr' :: VU.Vector Nuc -> DIM2 -> (Scalar (VU.Vector Nuc))
regionpr' inp (Z:.i:.j)
  = assert (let n = VU.length inp -1 in i<=j && i>=0 && j+1<=n+1)
  . Scalar $ VU.unsafeSlice i (j-i+1) inp
{-# INLINE regionpr' #-}

-- | Tests if (i,j) is a valid base pair. 

basepairing' inp (Z:.i:.j) = tf
  where p     = mkViennaPair (inp!(Z:.i), inp!(Z:.j-1))
        Z:.n  = snd . bounds $ inp
        tf    = i>=0 && j>0 && i<=j && j-1<=n && i<=n && p /= vpNS
{-# INLINE basepairing' #-}

stackpairing' inp k (Z:.i:.j) = tf
  where ps   = [ mkViennaPair (inp!(Z:.i+l), inp!(Z:.j-1-l)) | l <-[0..k-1] ]
        Z:.n = snd . bounds $ inp
        tf   = i>=k-1 && j>k-1 && i+k-1<=j-k+1 && j-k<=n && i+k-1<=n && all (/=vpNS) ps
{-# INLINE stackpairing' #-}

constrained cns ij = cns ij
{-# INLINE constrained #-}

-- |

reglen' :: Primary -> DIM2 -> Scalar Int
reglen' inp (Z:.i:.j)
  = assert (let (Z:.o,Z:.n) = bounds inp in i>=0 && (j<=n+1 || n== -1) && i<=j)
  . Scalar $ j-i
{-# INLINE reglen' #-}

-- |

reglenpl' :: Primary -> DIM2 -> Scalar (Nuc,Nuc,Int)
reglenpl' inp (Z:.i:.j)
  = assert (let (_,Z:.n) = bounds inp in i>0 && j<=n && i<=j)
  . Scalar $ (index inp (Z:.i-1),index inp (Z:.i),j-i)
{-# INLINE reglenpl' #-}

reglenpr' :: Primary -> DIM2 -> Scalar (Int,Nuc,Nuc)
reglenpr' inp (Z:.i:.j)
  = assert (let (_,Z:.n) = bounds inp in i>=0 && j<n && i<=j)
  . Scalar $ (j-i,index inp (Z:.j-1), index inp (Z:.j))
{-# INLINE reglenpr' #-}

-- | True, if the subword at ij is empty.

empty :: DIM2 -> Scalar Bool
empty (Z:.i:.j) = Scalar $ i==j

