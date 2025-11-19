# cubature

Cubature is an R package for adaptive multivariate integration over
hypercubes using deterministic and Monte Carlo methods. The package
provides wrappers around two C libraries:

- The cubature C library of [Steven G.
  Johnson](https://math.mit.edu/~stevenj/)
- The Cuba C library of [Thomas
  Hahn](https://wwwth.mpp.mpg.de/members/hahn/)

Both scalar and vectorized interfaces are available for all routines. R
users will gain a lot by using the vectorization as demonstrated in the
included vignette.

## Dedication

Version 2.0 is dedicated to the memory of Kiên Kiêu, one of the original
authors of `R2Cuba`. I (BN) never knew him, but this package derives
from his work on `R2Cuba`.
