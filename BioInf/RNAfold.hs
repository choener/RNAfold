{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RecordWildCards #-}

-- | RNAfold using the 2013 fusion framework.
--
-- TODO Split the interior loop calculations into different parts to speed up
-- calculations.

module BioInf.RNAfold where


import Biobase.Turner

import Debug.Trace




