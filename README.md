
ViennaRNA RNAfold v2, MFE variant
using the ADPfusion library



Introduction
============

This algorithm is the second, and much larger, test case for ADPfusion. We
implement "RNAfold v2" in the MFE variant using "-d2" dangles. Both a library
version and an executable are created. The "RNAFold" binary expects single
sequences, one per line. Backtracking tracks all co-optimal structures.



Installation
============

A simple "cabal update && cabal-dev install RNAFold" should be enough.



Runtime notes
=============

Using Haskell and ADPfusion, we come to within x3-x4 for this package. Between
the initial test case / submission (in 0.0.0.3) I have traded in some
performance improvements for much better readability in BioInf.RNAfold.Energy.
The C version of RNAfold employs some other methods to improve performance.
Consider:

base -~+ inner-1 +~- base
base -~+ inner-2 +~- base

where it is advantageous to calculate the outer basepair only once, not twice
as we are doing. It is probably better to try to improve the handling of
fusioned code and/or final assembler generation than finding calculations
common to different parts of CFG's.
