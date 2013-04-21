{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE RecordWildCards #-}

-- | Simple wrapper around the rnafold library.

module Main where



import System.Console.CmdArgs

import Biobase.Primary
import Biobase.Vienna
import BioInf.RNAfold
import qualified Biobase.Turner.Import as TI



data Options = Options
  { params :: String
  } deriving (Show,Data,Typeable)

options = Options
  { params = "./params" &= help "Turner 2004 RNA parameters (defaults to ./params)"
  }

main = do
  Options{..} <- cmdArgs options
  xs <- fmap lines getContents
  tm <- fmap turnerToVienna $ TI.fromDir params "" ".dat"
  mapM_ (run' tm) xs

run' tm inp = do
  print $ length inp
  print $ rnaFold tm (mkPrimary inp)

