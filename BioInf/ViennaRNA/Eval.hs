{-# LANGUAGE PatternGuards #-}

-- Direct evaluation of the energy of a given structure. The RNAfold-based
-- variant finds the optimal subset of base pairs that conform to the given
-- structure, this algorithm gives the energy of exactly the given structure.

module BioInf.ViennaRNA.Eval where

import Data.Vector.Fusion.Util (Id(..))
import qualified Data.Vector.Unboxed as VU
import Text.Printf

import Biobase.Primary
import Biobase.Secondary
import Biobase.Secondary.Diagrams
import Biobase.Vienna

import BioInf.ViennaRNA.Signature
import BioInf.ViennaRNA.Energy



rnaEval ener s d1s = flatten $ eval mfe ener s d1s

flatten :: SSTree PairIdx Deka -> (Deka, [String])
flatten = f where
  unDeka (Deka e) = e
  f (SSExt l e xs) =
    let etot = e + sum (map fst ys)
        ys   = map f xs
        here = printf "External loop: %d" (unDeka e)
    in  (etot, here : concatMap snd ys)
  f (SSTree p e xs) =
    let etot = e + sum (map fst ys)
        ys   = map f xs
    in  (etot, [])

eval :: Signature Id Deka Deka -> Vienna2004 -> Primary -> D1Secondary -> SSTree PairIdx Deka
eval efun ener s d1s = annotateWithEnergy t where
  t = d1sTree d1s
  (hairpin,interior,multi,blockStem,blockUnpair,compsBR,compsBC,structW,structCS,structWS,structOpen,h) = efun
  annotateWithEnergy :: SSTree PairIdx () -> SSTree PairIdx Deka
  annotateWithEnergy (SSExt l () xs) = SSExt l e (map annotateWithEnergy xs) where
    e = 0 -- TODO sum of all external loop energies
  annotateWithEnergy err@(SSTree (i,j) () xs)
    -- hairpin
    | null xs
    = SSTree (i,j) (hairpin ener si sii (VU.slice i (j-i) s) jjs sj) []
    -- interior loop
    | [SSTree (k,l) () _] <- xs
    = let kks = s VU.! (k-1); sll = s VU.! (l+1)
      in  SSTree (i,j) (interior ener si (VU.slice (i+1) (k-i-1) s) kks 0 sll (VU.slice (l+1) (j-l-1) s) sj) (map annotateWithEnergy xs)
    -- multibranched loop
    | otherwise
    = let e = multi ener si sii 0 0 jjs sj
            + sum (map bStem xs)
            + sum (map (\c -> blockUnpair ener c 0) cs)
          cs = [] -- TODO all unpaired nucleotides
          bStem (SSTree (k,l) () _) =
            let kks = s VU.! (k-1)
                sk  = s VU.! k
                sl  = s VU.! l
                sll = s VU.! (l+1)
            in  blockStem ener kks sk 0 sl sll
      in  SSTree (i,j) e (map annotateWithEnergy xs)
    where
      si  = s VU.! i
      sj  = s VU.! j
      sii = s VU.! (i+1)
      jjs = s VU.! (j-1)

