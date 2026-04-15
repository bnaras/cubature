# Build Notes

Personal notes for synchronizing to a new version of the underlying C
libraries.

## For libcubature (Steven G. Johnson's cubature)

Starting with R-package version 2.2.0, the upstream C library lives in
a separate repo
[`bnaras/libcubature`](https://github.com/bnaras/libcubature), forked
from [`stevengj/cubature`](https://github.com/stevengj/cubature) and
renamed to avoid colliding with the R package's own GitHub repo name.

The R package pulls it in as a submodule at `src/libcubature`,
pinned to the `R_pkg` branch. The `R_pkg` branch on
`bnaras/libcubature` carries these R-package-specific patches on top
of upstream `master`:

- `clencurt.h` regenerated with `M = 16` (not upstream's `M = 19`) to
  comply with CRAN package size limits
- `Makefile` replaced with a minimal one whose only target is
  `libcubature.a` (upstream's Makefile builds test harnesses instead)
- `_R_INTERFACE` preprocessor hook in `hcubature.c` so empty-heap
  errors go through `Rf_error()` when the library is compiled inside R
- Optional robust error estimation and new `NOT_CONVERGED` return
  code (new entry points `hcubature_robust` / `hcubature_v_robust` /
  `pcubature_robust` / `pcubature_v_robust`). See this package's
  `NEWS.md` and the "Robustness" section of the Get Started vignette
  for details. Upstream `stevengj/cubature` is dormant so we do not
  plan to upstream these changes.

Steps for updating to a new upstream release (should Steven ever cut
one, or for tracking upstream fixes):

1. In a clone of `bnaras/libcubature`: `git fetch origin` and merge
   (or rebase) upstream `master` into local `master`.
2. Rebase our `R_pkg` branch onto the new `master`. Expect merge
   conflicts in `hcubature.c`, `pcubature.c`, `cubature.h`, `Makefile`
   and `clencurt.h` — these are the R-package-specific patches.
3. Run `make libcubature.a` in `bnaras/libcubature` to verify the
   patched library still builds cleanly.
4. Force-push `R_pkg` back to `origin`.
5. In this R package clone, run
   `git submodule update --remote src/libcubature` to bump the
   submodule pin to the new `R_pkg` HEAD.
6. `R CMD INSTALL --preclean .` and `devtools::test()` to verify
   nothing broke.
7. Update `NEWS.md` to note the upstream-sync bump.

## For Cuba

Analogous arrangement: upstream lives at
[`bnaras/Cuba`](https://github.com/bnaras/Cuba), forked from Thomas
Hahn's tree so we can carry local patches. The R package uses it as a
submodule at `src/Cuba`, pinned to the `R_pkg` branch. Unix-like
systems are the default; Windows-specific fixes are isolated in the
`Cuba-win` directory and copied over via `Makefile.win`.

Update steps are the same shape as for libcubature: fetch upstream,
rebase `R_pkg`, push, then bump the submodule in this repo.

