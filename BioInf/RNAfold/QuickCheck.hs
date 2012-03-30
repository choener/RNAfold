
module BioInf.RNAfold.QuickCheck where

import Test.QuickCheck
import Test.QuickCheck.All

import ADP.Fusion.QuickCheck.Arbitrary -- hiding (options,customCheck,allProps)



-- * property checking

fCombined (i,j) = S.toList $ (,,) F.<<< fRegion #~~ fRegion ~~# fRegion F.... id $ Z:.i:.j

bCombined (i,j) = [ ( (i,k),(k,l),(l,j) )
                  | k <- [i..j]
                  , l <- [k..j]
                  , k-i >= 3
                  , j-l >= 3
                  , (k-i) + (j-l) >= 8
                  , (k-i) + (j-l) <= 32
                  ]

prop_Combined (Small i, Small j) = fCombined (i,j) == bCombined (i,j)

options = stdArgs {maxSuccess = 1000}

customCheck = quickCheckWithResult options

allProps = $forAllProperties customCheck

