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

import Debug.Trace



rnaEval ener s d1s = flatten $ eval mfe ener s d1s

flatten :: SSTree PairIdx Structure -> (Deka, [String])
flatten = f where
  unDeka (Deka e) = e
  f (SSExt l (External e) xs) =
    let etot = e + sum (map fst ys)
        ys   = map f xs
        here = printf "External loop: %d" (unDeka e)
    in  (etot, here : concatMap snd ys)
  f (SSTree _ (Hairpin  e l us r)          [] ) = (e, [printf "Hairpin loop: %d" (unDeka e)])
  f (SSTree _ (Interior e l ls ll rr rs r) [y]) =
    let etot = e + fst (f y)
    in  (etot, printf "Interior loop: %d" (unDeka e) : snd (f y))
  f (SSTree _ (Multi    e ll l r rr)       ys)  =
    let etot = e + sum (map (fst . f) ys)
    in  (etot, printf "Multi loop: %d %s" (unDeka e) (concatMap show [ll,l,r,rr]) : concatMap (snd . f) ys)
  {-
  f (SSTree p e xs) =
    let etot = e + sum (map fst ys)
        ys   = map f xs
        here
          | null xs   = printf "Hairpin loop: %d" (unDeka e)
          | [_] <- xs = printf "Interior loop: %d" (unDeka e)
          | otherwise = printf "Multibranched loop: %d" (unDeka e)
    in  (etot, here : concatMap snd ys)
    -}

data Structure
  = External Deka
  | Hairpin  Deka Nuc Primary Nuc
  | Interior Deka Nuc Primary Nuc Nuc Primary Nuc
  | Multi    Deka Nuc Nuc Nuc Nuc

eval :: Signature Id Deka Deka -> Vienna2004 -> Primary -> D1Secondary -> SSTree PairIdx Structure
eval efun ener s d1s = annotateWithEnergy t where
  t = d1sTree d1s
  (hairpin,interior,multi,blockStem,blockUnpair,compsBR,compsBC,structW,structCS,structWS,structOpen,h) = efun
  annotateWithEnergy :: SSTree PairIdx () -> SSTree PairIdx Structure
  annotateWithEnergy (SSExt l () xs) = SSExt l e (map annotateWithEnergy xs) where
    e = External 0 -- TODO sum of all external loop energies
  annotateWithEnergy err@(SSTree (i,j) () xs)
    -- hairpin
    | null xs
    = let pri = VU.slice (i+1) (j-i-1) s in SSTree (i,j) (Hairpin (hairpin ener si sii pri jjs sj) si pri sj) []
    -- interior loop
    | [SSTree (k,l) () _] <- xs
    = let kks = s VU.! k; sll = s VU.! l
          e   = interior ener si ls kks 0 sll rs sj
          ls  = VU.slice (i+1) (k-i-1) s
          rs  = VU.slice (l+1) (j-l-1) s
      in  SSTree (i,j) (Interior e si ls kks sll rs sj) (map annotateWithEnergy xs)
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
      in  SSTree (i,j) (Multi e si sii jjs sj) (map annotateWithEnergy xs)
    where
      si  = s VU.! i
      sj  = s VU.! j
      sii = s VU.! (i+1)
      jjs = s VU.! (j-1)

