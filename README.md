## cubature [![Build Status](](https://travis-ci.org/cubature) 

Cubature is an R package for adaptive multivariate integration over
hypercubes. The core of the package is the cubature C library
of [Steven G. Johnson](https://math.mit.edu/~stevenj/).

The package provides both `hcubature` and `pcubature` routines in
addition to a vector interface. R users will gain a lot by using the
vector interace and writing functions to take advantage of the vector
interface is quite easy. Refer to the vignette that provides timing
comparisons and also examples to emulate for using the vector
interface.
	
### Change log

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
