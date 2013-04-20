{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | RNAEval tool.

module Main where



import System.Console.CmdArgs

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
  mapM_ (run' tm) $ toPairs xs

toPairs (x1:x2:xs) = (x1,x2) : toPairs xs
toPairs [x] = error "single last line remaining"
toPairs [] = []

run' tm (inp,str) = do
  print $ length inp
  print $ rnaEval tm inp str



test inp str = do
  tm <- fmap turnerToVienna $ TI.fromDir "./params" "" ".dat"
  run' tm (inp,str)

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

