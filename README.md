
<!-- README.md is generated from the source: README.Rmd -->
cubature
========

[![Travis-CI Build Status](https://travis-ci.org/bnaras/cubature.svg?branch=master)](https://travis-ci.org/bnaras/cubature) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/cubature)](https://cran.r-project.org/package=cubature) [![Coverage Status](https://img.shields.io/codecov/c/github/bnaras/cubature/master.svg)](https://codecov.io/github/bnaras/cubature?branch=master)

Cubature is an R package for adaptive multivariate integration over
hypercubes. The core of the package is the cubature C library of
[Steven G. Johnson](https://math.mit.edu/~stevenj/).

The package provides both `hcubature` and `pcubature` routines in
addition to a vector interface.

R users will gain a lot by using the vector interace. Writing
functions to take advantage of the vector interface is quite
easy. Refer to the vignette that provides timing comparisons and
vectorization examples to emulate.
