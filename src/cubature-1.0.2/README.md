# cubature-1.0.2

This is a minimally modified copy of
[cubature-1.0.2 C library](http://ab-initio.mit.edu/wiki/index.php/Cubature)
needed for the for
[R package `cubature`](https://cran.r-project.org/package=cubature).

This C library has a large include file `clencurt.h` that violates
CRAN limitations on size. Therefore, I generated the file for M=16
instead of the default M=19. If you wish to use a higher M, you
have to link the clencurt_gen.c against the `long double` version
of the `fftw3` library. It would seem this would be easy to make
this customizable, but the R `fftw3` library is not distributed in a
linkable format. That is not the real impediment however; I have
not examined the effect of the `double` versus `long double`
in order to use it in production code. I might do so in the
future.

The specific modifications/additions I made are:

- This file

- Added a Makefile for the R package

- Added an `ifdef` macro to avoid the compilation warnings for a couple
  of lines of `hcubature.c` (line 920) and `pcubature.c` (line 310).

