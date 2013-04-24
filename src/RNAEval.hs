{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | RNAEval tool.

module Main where



import System.Console.CmdArgs

import Biobase.Primary
import Biobase.Secondary.Diagrams
import Biobase.Vienna
import qualified Biobase.Turner.Import as TI

import BioInf.ViennaRNA.Fold
import BioInf.ViennaRNA.Eval



data Options
  = Eval
      { params :: String
      }
  | ConstrainedFold
      { params :: String
      }
  deriving (Show,Data,Typeable)

oEval = Eval
  { params = "./params" &= help "Turner 2004 RNA parameters (defaults to ./params)"
  }

oConstrainedFold = ConstrainedFold
  { params = "../params"
  }

main = do
  o <- cmdArgs $ modes [oEval &= auto, oConstrainedFold]
  xs <- fmap lines getContents
  tm <- fmap turnerToVienna $ TI.fromDir (params o) "" ".dat"
  case o of
    Eval{..}            -> mapM_ (doEval tm) $ toPairs xs
    ConstrainedFold{..} -> mapM_ (doCF   tm) $ toPairs xs

toPairs (x1:x2:xs) = (x1,x2) : toPairs xs
toPairs [x] = error "single last line remaining"
toPairs [] = []

doEval tm (inp,str) = do
  print $ length inp
  print $ rnaEval tm (mkPrimary inp) (mkD1S str)

doCF tm (inp,str) = do
  print $ length inp
  print $ rnaFoldConstrained tm (mkPrimary inp) (mkD1S str)



test inp str = do
  tm <- fmap turnerToVienna $ TI.fromDir "./params" "" ".dat"
  doEval tm (inp,str)

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

