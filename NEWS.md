# cubature 2.1.1

- Fix documentation errors shown in CRAN checks
- Change default `absError` value in `hcubature`, `adaptIntegrate`,
  `pcubature` to `.Machine$double.eps * 10^2` to avoid nontermination,
  a breaking change unfortunately, but necessary
- Address issue raised by @jeroen about cross-compilation (thanks Jeroen).

# cubature 2.1.0

- Update cubature to 1.0.4 and Cuba to 4.2.2
- Added RANLIB call to `libcubature.a` ([Issue 40](https://github.com/bnaras/cubature/issues/40)). Thanks Sergei Fedorov.
- Simplified windows/non-windows organization, with default being windows

# cubature 2.0.4.6

- Fix CRAN C23 issues using system requirements USE_C17
- Bug fix: increment correct count of number of evaluations. Thanks Jan Meis.
- Removed all solaris and unneeded Rcpp flags
- Copied correct includes to inst/include for linking

# cubature 2.0.4.5

- Fix math typesetting in documentation and use roxygen markdown tags wherever
  possible.

# cubature 2.0.4.4

- Fix unprotected arguments ([Issue 34](https://github.com/bnaras/cubature/issues/34)). Thanks Manuel Koller.
- Fix documentation ([Issue 32](https://github.com/bnaras/cubature/issues/32)). Thanks Emanuele Guidotti.

# cubature 2.0.4.3

- Unreleased

# cubature 2.0.4.2

- Fix `Makevars` to refer to `CC` not `C` ([Issue 33](https://github.com/bnaras/cubature/issues/33))
- Replace `src/Cuba-4.2-win/common/Random.c:105` with equivalent to ensure gcc 8.x erroneous warning is avoided on Windows
- Add `rmarkdown` to suggests in description

# cubature 2.0.4.1

- Pass `$(AR)`, `$(ARFLAGS)`, `$(RANLIB)` to `make` so that LTO checks
  pick up correct plugins. Also add `cleanup` script.  (Thanks,
  Prof. Brian Ripley)

# cubature 2.0.4

- `cubintegrate` now matches method via `match.arg` ([Issue
  25](https://github.com/bnaras/cubature/issues/25))

- Address `gcc` version 10.0 changes due to `-fno-common` default
  setting.

# cubature 2.0.3

- Fixed up stack overrun in `Cuba-4.2/src/divonne/Split.c` (lines
  119--128 utilizing a flag for first time through loop
- More cleanup of `Makevars`

# cubature 2.0.2

- Fixed up uninitialized count for `hcubature` (Thanks, Ehsan
  Masoudi)
- Tentative fix of infinite value for `minfluct` in
  `Cuba-4.2/src/suave/Integrate.c` in lines 197--205. (Thanks,
  Prof. Brian Ripley) 
- Miscellaneous fixes for compilation on various platforms.
  * Cleaned up `Makevars` and `Makevars.win` to remove some unused
    flags. Tentative fix for solaris in `Makevars` and
    `Cuba-4.2/src/common/stddecl.c`.
  * Removed printing of `statefile` in all the `src/*/Integrate.c`
    routines to avoid segfault on solaris.
  * Added `ret_code` for return value in call to `getloadavg` in
    `Cuba-4.2/src/common/Fork.c`. Also added an `#ifdef SOLARIS` for
    including appropriate header under solaris.
  * Fixed up call to `MASTER` variadic macro by including dummy argument
    0 in `Cuba-4.2/src/common/Fork.c`, line 146.
  * Replaced multi-statement macro with a `do { } while(0)` hack.
  * Reworked code around `Cuba-4.2/src/divonne/Rule.c` checks for
    key values (lines 593 to 600).
  * Fixed up call to `WORKER` variadic macro by including dummy argument
    0 in `Cuba-4.2/src/common/Parallel.c`, line 393.
  * Moved `#ifdef FRAMECOPY` outside in
    `Cuba-4.2/src/common/Parallel.c` to avoid embedding directives in
    macro args, lines 70-80.

# cubature 2.0

- Major update. Integrates [Cuba 4.2](https://feynarts.de/cuba/)
  library
- Allows for finite or infinite limits of integration (courtesy of
  Simon Gaure)
- Deprecates `doChecking` argument. Now it does nothing and will be
  removed in future versions.

# cubature 1.4-1

- Fixed up `scale` argument to tests to conform to pass checks.
- Removed any reference to orphaned package `R2Cuba` and updated
  vignette with information on imminent 2.0 release.

# cubature 1.4

- Fixed up the C call so that it is re-entrant (brought to my
  attention by Pierre de Villemereuil). This should be considered a
  bug fix!

- Corrected private notes

# cubature 1.3-13

- Generated package docs using pkgdown
- Synced up to cubature-1.0.3 
- Added LICENSE file ([request of Nick Youngblut](https://github.com/bnaras/cubature/issues/12))

# cubature 1.3-12

- Minor typographical and documentation fixes

# cubature 1.3-11

- Merged Manuel Koller's registration for C code and vignette fix
  for NA.
- Moved cubature header and exp_cubature headers to `inst/include`
  for linking to other packages

# cubature 1.3-10

- Renamed `Readme.md` to `README.md`
- Removed references to ab-initio website that caused some hassle
  due to misconfigured site.

# cubature 1.3-9

# cubature 1.3-8

# cubature 1.3-7

- Registered .Call stuff and removed microbenchmark suggestion in
  favor of benchr.

# cubature 1.3-4

- Generated smaller pcubature header (`clencurt.h`) using
  M = 16 and put back C cubature source in tree.

# cubature 1.3-3

- Moved the cubature-1.0.2 C library to github to avoid
  hitting the CRAN size limit.

# cubature 1.3-2

- Added vignette showing huge gains due to vectorization

# cubature 1.3-1

- Added pcubature
- Added vector versions

# cubature 1.3-0

- Moved to SGJ cubature-1.0.2 version.
- `adaptIntegrate` and `hcubature` are aliases
- `hcubature` function gains a norm argument that is set to a
  sensible default so as not to affect depending packages

# cubature 1.2-0

- Moved to `Rcpp` framework
- Added tests in preparation for SGJ cubature-1.0.2 version

# cubature 1.1-3

- Roxygenized in preparation for upgrade to newer version of
  cubature library on abinitio website

# cubature 1.1-2

- Registered native cubature functions `adapt_integrate` and
  `adapt_integrate_v` so that they are directly callable from C
  (courtesy of Simen Gaure)

# cubature 1.1-2

- Fixed typo in doc for function `adaptIntegrate`; default value for
  `doChecking` was incorrectly stated as `TRUE`

# cubature 1.1-1

- Added `doChecking` argument (default `FALSE`) to save some
  computation time in evaluating integrand (9% speedup).

# cubature 1.1

- Synced up to SGJ cubature routines dated 2010-10-18 on [his
  website](https://github.com/stevengj/cubature)
- Bugfix: potential memory leak fixed up in heap routine (my oversight!)
- Routine `adaptIntegrate` gains ... argument (request of Baptiste Auguie)
- Corrected radius constant in `testFn2` to match cubature output exactly

# cubature 1.0
  - Original version of package based on Steven G. Johnson's cubature
	routines at http://ab-initio.mit.edu/wiki/index.php/Cubature

