## Resubmission

This package was archived on 2022-06-14 due to C compilation failures
caused by deprecated R API macros (`Calloc`/`Free`). All issues have
been resolved:

- Replaced `Calloc`/`Free` with `R_Calloc`/`R_Free` throughout
- Added native routine registration (`src/init.c`)
- Modernized NAMESPACE with explicit `importFrom` declarations
- Updated DESCRIPTION to `Authors@R` format

Additionally, this version adds:

- OpenMP parallelization for `learnPattern` and `computeSimilarity`
  (with `nthreads` parameter, graceful fallback when OpenMP unavailable)
- New `discriminativePatterns()` function for class-discriminative
  pattern identification and visualization
- Improved `visualizePattern` plot aesthetics

## Test environments

- Ubuntu 22.04 LTS, R 4.5.2, GCC 11.4.0 (with OpenMP)

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: New submission / Package was archived on CRAN
  - This is a resubmission with all prior issues fixed.
