
module BioInf.ViennaRNA.Signature where

import Data.Vector.Fusion.Stream.Monadic as SM

import Biobase.Primary
import Biobase.Vienna



type Signature m a r =
  -- weak / hairpin
  ( Vienna2004 -> Nuc -> Nuc -> Primary -> Nuc -> Nuc -> a
  -- weak / interior
  , Vienna2004 -> Nuc -> Primary -> Nuc -> a -> Nuc -> Primary -> Nuc -> a
  -- weak / multibranch
  , Vienna2004 -> Nuc -> Nuc -> a -> a -> Nuc -> Nuc -> a
  -- block / multistem
  , Vienna2004 -> Nuc -> Nuc -> a -> Nuc -> Nuc -> a
  -- block / unpaired
  , Vienna2004 -> Nuc -> a -> a
  -- comps / block region
  , Vienna2004 -> a -> Primary -> a
  -- comps / block comps
  , Vienna2004 -> a -> a -> a
  -- struct / weak
  , Vienna2004 -> a -> a
  -- struct / char-struct
  , Vienna2004 -> Nuc -> a -> a
  -- struct / weak-struct
  , Vienna2004 -> a -> a -> a
  -- struct / open
  , Vienna2004 -> Primary -> a
  -- all / objective
  , Stream m a -> m r
  )


