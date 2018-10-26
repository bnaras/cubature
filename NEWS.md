# cubature 2.0

- Major update. Integrates [Cuba 4.2](http://www.feynarts.de/cuba/)
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

