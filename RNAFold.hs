
-- | Simple wrapper around STrnafold

module Main where

import BioInf.RNAfold



main = do
  xs <- fmap lines getContents
  mapM_ doRNAfold xs

doRNAfold inp = do
  let (e,bs) = testRNAfold inp
  putStrLn inp
  print e
  mapM_ putStrLn bs
