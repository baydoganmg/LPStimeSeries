# LPStimeSeries TODO

## Package Revival (Completed)
- [x] Fix C compilation: `Calloc`/`Free` -> `R_Calloc`/`R_Free`
- [x] Add native routine registration (`src/init.c`)
- [x] Add missing `importFrom` declarations in NAMESPACE
- [x] Update DESCRIPTION metadata (`Authors@R`, R >= 3.5.0, etc.)
- [x] Modernize CITATION to `bibentry` format
- [x] Replace dead homepage URLs with DOI in man pages
- [x] Remove stale MD5 file
- [x] Pass `R CMD check` with Status: OK on R 4.5.2

## Improvements (from inst/NEWS TODO)

### 1. Speed up `tunelearnPattern`
**Status:** Investigated, deferred.
**Finding:** Benchmarked on GunPoint — training (`learnPattern`) accounts for ~68% of time, similarity computation only ~22%. Switching to representation-based similarity saves ~4% of total time. Not worth the added complexity for now. Revisit if tuning becomes a bottleneck on larger datasets.

### 2. Simplify tree size control
Currently grows full tree by `nodesize` then truncates by `maxdepth`, which is computationally wasteful. Plan was to drop depth parameter and control tree size by a single parameter.

### 3. Parallelization
**Status:** Completed via OpenMP.
Added `nthreads` parameter to `learnPattern` and `computeSimilarity`. Tree building and similarity computation are parallelized with per-thread RNG (xoshiro256**). Benchmarks show near-linear scaling (e.g., ~3.4x at 4 threads, ~5.9x at 8 threads on 500×150 dataset with 500 trees). Falls back gracefully to single-threaded when OpenMP is unavailable.

### 4. Fair segment sampling
Observations near the start and end of time series are underrepresented because segment start positions are sampled uniformly from `[0, length - segment_length)`. A fair sampling strategy is needed.

### 5. Better pattern visualization
**Status:** Completed.
- Improved `visualizePattern`: full series as faint background line, shaded predictor/target time windows, connected line segments with points instead of raw scatter.
- Added `discriminativePatterns` function: ranks terminal nodes by chi-squared-like class discriminability score, plots top patterns with per-class colored background series, shaded windows, per-series and mean pattern lines.
- Both functions handle non-contiguous time segments correctly (no lines drawn across gaps).

### 6. Variable-length time series support
Currently requires UCR format (equal-length rows). LPS can handle variable-length series but the input format needs to change.

### 7. Multivariate time series support
LPS can work for multivariate time series similarity. Current version needs modifications.

### 8. Categorical time series support
LPS can work for categorical time series (e.g., DNA sequences). Current version requires modification for categorical variables.
