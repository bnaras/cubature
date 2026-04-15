# Critical Review: cubature R Package

**Date**: 2026-01-28 **Status**: Most issues addressed in v2.2.0

------------------------------------------------------------------------

## Strengths

1.  **Solid External Library Integration**: Clean wrapping of two
    well-established C libraries (cubature and Cuba) via Rcpp
2.  **Multi-platform Support**: Build configurations for Unix, Windows,
    and WebAssembly
3.  **Good CI/CD**: GitHub Actions for R-CMD-check across multiple R
    versions
4.  **Comprehensive Test Coverage**: Tests for both scalar and
    vectorized interfaces
5.  **Detailed Changelog**: NEWS.md shows careful attention to backwards
    compatibility and bug fixes

------------------------------------------------------------------------

## Issues and Status

### 1. Build System Complexity

- ~~**Timestamp-based builds**: Using `.ts` files for incremental builds
  is fragile. If a source file in libcubature/ or Cuba/ changes, the
  timestamp won’t update unless you manually delete the `.ts` file.~~
  - **Status**: FIXED - Added source file dependencies to all Makevars
    files
- **Cuba configure step**: Running `./configure` at build time
  (`Makevars:21`) can fail on restricted build environments and adds
  unpredictability.
  - **Status**: Not addressed (would require significant restructuring)
- **`SystemRequirements: GNU make and USE_C17`**: The USE_C17
  requirement suggests workarounds for C23 issues.
  - **Status**: Not addressed (required for CRAN compatibility)

### 2. Code Quality Issues

- ~~**Deprecated parameter still present** (`hcubature.R`): `doChecking`
  parameter documented as ignored since 2.0~~
  - **Status**: FIXED - Removed in v2.2.0
- ~~**Typo in default_args()** (`cubintegrate.R:60`): `sauve` instead of
  `suave`~~
  - **Status**: FIXED in v2.2.0
- ~~**Inconsistent default absError** between `hcubature` and
  `pcubature`~~
  - **Status**: FIXED - Both now use `.Machine$double.eps * 10^2`
- **Hard-coded PACKAGE argument**: Using `PACKAGE="cubature"` in
  [`.Call()`](https://rdrr.io/r/base/CallExternal.html) vs generated
  `RcppExports.R`
  - **Status**: Not addressed (minor, works correctly)

### 3. Testing Gaps

- ~~**No negative tests**: No tests for invalid inputs, error handling,
  edge cases~~
  - **Status**: FIXED - Added comprehensive test suite:
    - `test_input_validation.R` - Invalid input rejection
    - `test_error_handling.R` - Error propagation
    - `test_boundaries.R` - Edge cases (maxEval, infinite limits)
    - `test_consistency.R` - Method agreement
    - `test_regression.R` - Known reference values
- **Commented-out tests** (`test_cuba.R`): Some tests commented out as
  flaky
  - **Status**: Not addressed (may indicate underlying issues)

### 4. Documentation Issues

- ~~**Typo in return value docs** (`cubintegrate.R:125`): `forcCuba`
  should be `for Cuba`~~
  - **Status**: FIXED in v2.2.0
- ~~**Outdated example URLs**: References to ab-initio.mit.edu~~
  - **Status**: FIXED - Updated to GitHub URL
- ~~**Examples wrapped in `\dontrun{}`**: All examples skipped during R
  CMD check~~
  - **Status**: FIXED - Added quick examples that run, converted rest to
    `\donttest{}`

### 5. API Inconsistencies

- **Parameter naming**: `hcubature` uses `lowerLimit/upperLimit` while
  `cubintegrate` uses `lower/upper`
  - **Status**: Not addressed (would break backwards compatibility)
- **Return value naming**: `functionEvaluations` vs `neval`
  - **Status**: Not addressed (would break backwards compatibility)
- **fDim vs nComp**: Different names for same concept
  - **Status**: Not addressed (would break backwards compatibility)

### 6. Potential Memory/Safety Issues

- **R callback from C**: No visible error handling for R function errors
  in `fWrapper`
  - **Status**: Not addressed (documented as known limitation in tests)
- **No input validation in C++**: Rcpp wrappers trust R-side validation
  - **Status**: Not addressed (documented as known limitation)

### 7. Git Submodule Concerns

- **Cuba submodule points to fork**: May diverge from upstream
  - **Status**: Not addressed (intentional for R compatibility)
- ~~**cubature-1.0.4 is not a submodule**: Direct copy with
  modifications~~
  - **Status**: FIXED in v2.2.0 — `src/libcubature/` is now a git
    submodule pointing at `bnaras/libcubature@R_pkg`, paralleling the
    `src/Cuba` → `bnaras/Cuba@R_pkg` arrangement
  - **Status**: Not addressed (documented in NOTES.md)

### 8. CI/CD

- ~~**Missing Windows CI**: Only macOS and Ubuntu tested~~
  - **Status**: FIXED - Added Windows to CI matrix

------------------------------------------------------------------------

## Summary of Changes in v2.2.0

1.  Fixed typo `sauve` → `suave` in
    [`default_args()`](reference/default_args.md)
2.  Fixed typo `forcCuba` → `for Cuba` in documentation
3.  Harmonized `absError` defaults between hcubature and pcubature
4.  Updated outdated URL in examples
5.  Removed deprecated `doChecking` parameter (breaking change)
6.  Added Windows to CI test matrix
7.  Converted examples from `\dontrun{}` to `\donttest{}` with quick
    runnable examples
8.  Fixed timestamp-based builds to properly track source file
    dependencies
9.  Added comprehensive test suite (5 new test files, 78 new tests)
10. Regenerated documentation

------------------------------------------------------------------------

## Remaining Items (Not Addressed)

These items were not addressed due to backwards compatibility concerns,
low priority, or being inherent to the design:

- Cuba configure step at build time (would require significant
  restructuring)
- API naming inconsistencies (would break existing code)
- C++ layer error handling (documented limitation)
- Commented-out flaky tests (may need investigation)
