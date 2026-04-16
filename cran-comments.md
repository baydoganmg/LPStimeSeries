## Resubmission (2026-04-17)

In response to Benjamin Altmann's review, this version:

- Adds a methods reference to the DESCRIPTION field in the required
  form: Baydogan and Runger (2016) <doi:10.1007/s10618-015-0425-y>.
- Expands `Authors@R` to credit all contributors and copyright holders
  of adapted code. The C sources in `src/regTree.c`, `src/regrf.c`,
  `src/rfutils.c`, and `src/rf.h` are adapted from the 'randomForest'
  package (R port by Andy Liaw and Matthew Wiener, based on original
  Fortran code by Leo Breiman and Adele Cutler; copyright held by
  Merck & Co., Inc.). All four individuals are now listed with `ctb`
  role, and Merck & Co., Inc. with `cph` role. The original
  copyright and GPL notices remain intact at the top of each C file.

## Previous submission notes

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
