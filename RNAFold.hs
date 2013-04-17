
-- | Simple wrapper around STrnafold

module Main where

import BioInf.RNAfold
import qualified Biobase.Turner.Import as TI
import Biobase.Vienna



main = do
  xs <- fmap lines getContents
  tm <- fmap turnerToVienna $ TI.fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
  mapM_ (run' tm) xs

run' tm inp = do
  print $ length inp
  print $ rnaFold tm inp

{-
doRNAfold inp = do
  let (e,bs) = testRNAfold inp
  putStrLn inp
  print e
  mapM_ putStrLn bs
-}
