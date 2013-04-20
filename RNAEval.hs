
-- | Simple wrapper around STrnafold

module Main where

import BioInf.RNAfold
import qualified Biobase.Turner.Import as TI
import Biobase.Vienna



main = do
  xs <- fmap lines getContents
  tm <- fmap turnerToVienna $ TI.fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
  mapM_ (run' tm) $ toPairs xs

test inp str = do
  tm <- fmap turnerToVienna $ TI.fromDir "/home/choener/Documents/Workdata/TurnerRNA2004/RNA" "" ".dat"
  run' tm (inp,str)

toPairs (x1:x2:xs) = (x1,x2) : toPairs xs
toPairs [x] = error "single last line remaining"
toPairs [] = []

run' tm (inp,str) = do
  print $ length inp
  print $ rnaEval tm inp str

tests = mapM_ (uncurry test)
  [ ( "CCUGACUGGCGUUGACAUAUGGUU"
    , ".......(((((......)).)))"
    )
  , ( "CUGGGGGUGACAUCCCCCC"
    , "..(((((......)).)))"
    )
  , ( "GGCGUUGACAUAUGGUU"
    , "(((((......)).)))"
    )
  , ( "GGGGUUGACAUACCCCC"
    , "(((((......)).)))"
    )
  , ( "GGCGUUGACAUAUGUU"
    , "(((((......)))))"
    )
  , ( "GGGGGUGACAUCCCCC"
    , "(((((......)))))"
    )
  , ( "GGGGGUGACCCCC"
    , "(((((...)))))"
    )
  , ( "CCUGACUGGCGUUGACAUAUGGUUGCUUGAGCGUAGCCAGGUGUUGGUGGUCCAGUGCAUCAAGGUGCCGUCGGAUCGGAUACUUGGCUUUGCUUAGAUU"
    , ".......(((((......)).)))(.(((((((.(((((((.(((((((.....((((......)))))))).)))......))))))).))))))).)."
    )
  ]


{-
doRNAfold inp = do
  let (e,bs) = testRNAfold inp
  putStrLn inp
  print e
  mapM_ putStrLn bs
-}
