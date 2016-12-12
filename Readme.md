cubature
========

Updated this package in several ways:

- Now using  `Rcpp` (as of version 1.2)
- Synced up to newer version of
  [cubature](http://ab-initio.mit.edu/wiki/index.php/Cubature) library
  (version 1.0.2) by Steven Johnson so that both `hcubature`
  `pcubature` are now available.
- Added vectorized interface to `hcubature` and `pcubature`
- Added a vignette to show the great gains from vectorization
- Added tests to avoid screw-ups
- Fixed up clencurt.h so that CRAN size limitations are not violated
