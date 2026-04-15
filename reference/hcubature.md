# Adaptive multivariate integration over hypercubes (hcubature and pcubature)

The function performs adaptive multidimensional integration (cubature)
of (possibly) vector-valued integrands over hypercubes. The function
includes a vector interface where the integrand may be evaluated at
several hundred points in a single call.

## Usage

``` r
hcubature(
  f,
  lowerLimit,
  upperLimit,
  ...,
  tol = 1e-05,
  fDim = 1,
  maxEval = 0,
  absError = .Machine$double.eps * 10^2,
  vectorInterface = FALSE,
  norm = c("INDIVIDUAL", "PAIRED", "L2", "L1", "LINF"),
  robust = FALSE
)

pcubature(
  f,
  lowerLimit,
  upperLimit,
  ...,
  tol = 1e-05,
  fDim = 1,
  maxEval = 0,
  absError = .Machine$double.eps * 10^2,
  vectorInterface = FALSE,
  norm = c("INDIVIDUAL", "PAIRED", "L2", "L1", "LINF"),
  robust = FALSE
)
```

## Arguments

- f:

  The function (integrand) to be integrated

- lowerLimit:

  The lower limit of integration, a vector for hypercubes

- upperLimit:

  The upper limit of integration, a vector for hypercubes

- ...:

  All other arguments passed to the function f

- tol:

  The maximum tolerance, default 1e-5.

- fDim:

  The dimension of the integrand, default 1, bears no relation to the
  dimension of the hypercube

- maxEval:

  The maximum number of function evaluations needed, default 0 implying
  no limit. Note that the actual number of function evaluations
  performed is only approximately guaranteed not to exceed this number.

- absError:

  The maximum absolute error tolerated, default
  `.Machine$double.eps * 10^2`.

- vectorInterface:

  A flag that indicates whether to use the vector interface and is by
  default FALSE. See details below

- norm:

  For vector-valued integrands, `norm` specifies the norm that is used
  to measure the error and determine convergence properties. See below.

- robust:

  Logical; defaults to `FALSE`. If `TRUE`, `hcubature` uses a more
  conservative error-estimation strategy combining a degree-3 null rule
  with the existing degree-5 null rule via a Cuhre-style decay-check
  formula, a parent-child consistency check at the subdivision step, and
  a peak-density suspect-region detector. The robust path is opt-in
  because it can increase the number of function evaluations (by 0% on
  most smooth integrands and up to ~30% on sharply-peaked ones), and it
  addresses silent "fool's-convergence" failures observed on non-smooth
  integrands with thin or localized support (e.g. a `(max(L(x), 0))^2`
  kink, a nearly-step logistic, or a very-narrow Gaussian inside an
  oversized bounding box). It is not a complete cure for all such cases
  — if the integrand's support is much smaller than the inter-sample
  spacing of the base rule, no amount of error-model cleverness can
  detect it; for those cases shrink the integration domain or use
  `cubintegrate(method = "cuhre")`. Only applies to `hcubature`;
  `pcubature` ignores it.

## Value

The returned value is a list of four items:

- integral:

  the value of the integral (a numeric vector of length `fDim` for
  vector integrands)

- error:

  the estimated absolute error on `integral`, as reported by the
  underlying routine. For single-integrand problems this is a single
  number; for vector integrands it is one error estimate per component.
  Callers who supplied a tolerance (`tol` for relative or `absError` for
  absolute) should treat the integral as reliable only when `error`
  satisfies that tolerance.

- functionEvaluations:

  the number of times the integrand was evaluated. Use this to budget
  `maxEval` across calls.

- returnCode:

  an integer status code from the underlying C routine, taking one of
  three values:

  `0` (success)

  :   the integrator reached the requested tolerance. The result can be
      trusted to `error` precision.

  `1` (failure)

  :   a hard internal error — allocation failure, invalid parameters, an
      integrand that returned a non-zero error code, or similar. In this
      case `integral` and `error` should be ignored; the call has not
      produced a meaningful result.

  `2` (not converged)

  :   the `maxEval` budget was exhausted before the requested tolerance
      was met. The `integral` and `error` fields still contain the best
      available estimate and its reported error bound, but the result
      has not been verified to satisfy the caller's tolerance and should
      be treated as best-effort. The R-level wrappers also emit a
      [`warning()`](https://rdrr.io/r/base/warning.html) in this case;
      wrap the call in
      [`suppressWarnings()`](https://rdrr.io/r/base/warning.html) to
      silence the warning if the best-effort value is acceptable for
      your use. This return code is new in cubature 2.2.0; prior
      versions silently reported `0` (success) on budget exhaustion,
      with the caller expected to compare `error` against their own
      tolerance. The symbolic constants `CUBATURE_SUCCESS`,
      `CUBATURE_FAILURE`, and `CUBATURE_NOT_CONVERGED` are exported from
      the upstream C header for compiled clients.

## Details

The function merely calls Johnson's C code and returns the results.

One can specify a maximum number of function evaluations (default is 0
for no limit). Otherwise, the integration stops when the estimated error
is less than the absolute error requested, or when the estimated error
is less than tol times the integral, in absolute value, or the maximum
number of iterations is reached (see parameter info below), whichever is
earlier.

For compatibility with earlier versions, the `adaptIntegrate` function
is an alias for the underlying `hcubature` function which uses
h-adaptive integration. Otherwise, the calling conventions are the same.

We highly recommend referring to the vignette to achieve the best
results!

The `hcubature` function is the h-adaptive version that recursively
partitions the integration domain into smaller subdomains, applying the
same integration rule to each, until convergence is achieved.

The p-adaptive version, `pcubature`, repeatedly doubles the degree of
the quadrature rules until convergence is achieved, and is based on a
tensor product of Clenshaw-Curtis quadrature rules. This algorithm is
often superior to h-adaptive integration for smooth integrands in a few
(\<=3) dimensions, but is a poor choice in higher dimensions or for
non-smooth integrands. Compare with `hcubature` which also takes the
same arguments.

The vector interface requires the integrand to take a matrix as its
argument. The return value should also be a matrix. The number of points
at which the integrand may be evaluated is not under user control: the
integration routine takes care of that and this number may run to
several hundreds. We strongly advise vectorization; see vignette.

The `norm` argument is irrelevant for scalar integrands and is ignored.
Given vectors \\v\\ and \\e\\ of estimated integrals and errors therein,
respectively, the `norm` argument takes on one of the following values:

- `INDIVIDUAL`:

  Convergence is achieved only when each integrand (each component of
  \\v\\ and \\e\\) individually satisfies the requested error tolerances

- `L1`, `L2`, `LINF`:

  The absolute error is measured as \\\|e\|\\ and the relative error as
  \\\|e\|/\|v\|\\, where \\\|...\|\\ is the \\L_1\\, \\L_2\\, or
  \\L\_{\infty}\\ norm, respectively

- `PAIRED`:

  Like `INDIVIDUAL`, except that the integrands are grouped into
  consecutive pairs, with the error tolerance applied in an \\L_2\\
  sense to each pair. This option is mainly useful for integrating
  vectors of complex numbers, where each consecutive pair of real
  integrands is the real and imaginary parts of a single complex
  integrand, and the concern is only the error in the complex plane
  rather than the error in the real and imaginary parts separately

## Author

Balasubramanian Narasimhan

## Examples

``` r
## Simple example: integrate product of cosines over [0,1]^2
testFn0 <- function(x) prod(cos(x))
hcubature(testFn0, rep(0, 2), rep(1, 2), tol = 1e-4)
#> $integral
#> [1] 0.7080734
#> 
#> $error
#> [1] 1.709434e-05
#> 
#> $functionEvaluations
#> [1] 17
#> 
#> $returnCode
#> [1] 0
#> 

## Simple polynomial: integral should be exactly 1
testFn3 <- function(x) prod(2 * x)
hcubature(testFn3, rep(0, 3), rep(1, 3), tol = 1e-4)
#> $integral
#> [1] 1
#> 
#> $error
#> [1] 2.220446e-16
#> 
#> $functionEvaluations
#> [1] 33
#> 
#> $returnCode
#> [1] 0
#> 

# \donttest{
## More examples (not run during R CMD check)

## Product of cosines with pcubature
pcubature(testFn0, rep(0, 2), rep(1, 2), tol = 1e-4)
#> $integral
#> [1] 0.7080734
#> 
#> $error
#> [1] 1.518315e-07
#> 
#> $functionEvaluations
#> [1] 81
#> 
#> $returnCode
#> [1] 0
#> 

M_2_SQRTPI <- 2/sqrt(pi)

## Test function 1
## Compare with original cubature result of
## ./cubature_test 3 1e-4 1 0
## 3-dim integral, tolerance = 0.0001
## integrand 1: integral = 1.00001, est err = 9.67798e-05, true err = 9.76919e-06
## #evals = 5115

testFn1 <- function(x) {
  val <- sum (((1-x) / x)^2)
  scale <- prod(M_2_SQRTPI/x^2)
  exp(-val) * scale
}

hcubature(testFn1, rep(0, 3), rep(1, 3), tol=1e-4)
#> $integral
#> [1] 1.00001
#> 
#> $error
#> [1] 9.677977e-05
#> 
#> $functionEvaluations
#> [1] 5115
#> 
#> $returnCode
#> [1] 0
#> 
pcubature(testFn1, rep(0, 3), rep(1, 3), tol=1e-4)
#> $integral
#> [1] NaN
#> 
#> $error
#> [1] 0
#> 
#> $functionEvaluations
#> [1] 27
#> 
#> $returnCode
#> [1] 0
#> 

##
## Test function 2
## Compare with original cubature result of
## ./cubature_test 2 1e-4 2 0
## 2-dim integral, tolerance = 0.0001
## integrand 2: integral = 0.19728, est err = 1.97261e-05, true err = 4.58316e-05
## #evals = 166141

testFn2 <- function(x) {
  ## discontinuous objective: volume of hypersphere
  radius <- as.double(0.50124145262344534123412)
  ifelse(sum(x*x) < radius*radius, 1, 0)
}

hcubature(testFn2, rep(0, 2), rep(1, 2), tol=1e-4)
#> $integral
#> [1] 0.1972834
#> 
#> $error
#> [1] 1.972631e-05
#> 
#> $functionEvaluations
#> [1] 164543
#> 
#> $returnCode
#> [1] 0
#> 
pcubature(testFn2, rep(0, 2), rep(1, 2), tol=1e-4)
#> $integral
#> [1] 0.1973732
#> 
#> $error
#> [1] 8.220791e-06
#> 
#> $functionEvaluations
#> [1] 1052929
#> 
#> $returnCode
#> [1] 0
#> 

##
## Test function 3
## Compare with original cubature result of
## ./cubature_test 3 1e-4 3 0
## 3-dim integral, tolerance = 0.0001
## integrand 3: integral = 1, est err = 0, true err = 2.22045e-16
## #evals = 33

testFn3 <- function(x) {
  prod(2*x)
}

hcubature(testFn3, rep(0,3), rep(1,3), tol=1e-4)
#> $integral
#> [1] 1
#> 
#> $error
#> [1] 2.220446e-16
#> 
#> $functionEvaluations
#> [1] 33
#> 
#> $returnCode
#> [1] 0
#> 
pcubature(testFn3, rep(0,3), rep(1,3), tol=1e-4)
#> $integral
#> [1] 1
#> 
#> $error
#> [1] 0
#> 
#> $functionEvaluations
#> [1] 27
#> 
#> $returnCode
#> [1] 0
#> 

##
## Test function 4 (Gaussian centered at 1/2)
## Compare with original cubature result of
## ./cubature_test 2 1e-4 4 0
## 2-dim integral, tolerance = 0.0001
## integrand 4: integral = 1, est err = 9.84399e-05, true err = 2.78894e-06
## #evals = 1853

testFn4 <- function(x) {
  a <- 0.1
  s <- sum((x - 0.5)^2)
  (M_2_SQRTPI / (2. * a))^length(x) * exp (-s / (a * a))
}

hcubature(testFn4, rep(0,2), rep(1,2), tol=1e-4)
#> $integral
#> [1] 1.000003
#> 
#> $error
#> [1] 9.843987e-05
#> 
#> $functionEvaluations
#> [1] 1853
#> 
#> $returnCode
#> [1] 0
#> 
pcubature(testFn4, rep(0,2), rep(1,2), tol=1e-4)
#> $integral
#> [1] 1
#> 
#> $error
#> [1] 4.415857e-09
#> 
#> $functionEvaluations
#> [1] 4225
#> 
#> $returnCode
#> [1] 0
#> 

##
## Test function 5 (double Gaussian)
## Compare with original cubature result of
## ./cubature_test 3 1e-4 5 0
## 3-dim integral, tolerance = 0.0001
## integrand 5: integral = 0.999994, est err = 9.98015e-05, true err = 6.33407e-06
## #evals = 59631

testFn5 <- function(x) {
  a <- 0.1
  s1 <- sum((x - 1/3)^2)
  s2 <- sum((x - 2/3)^2)
  0.5 * (M_2_SQRTPI / (2. * a))^length(x) * (exp(-s1 / (a * a)) + exp(-s2 / (a * a)))
}

hcubature(testFn5, rep(0,3), rep(1,3), tol=1e-4)
#> $integral
#> [1] 0.9999937
#> 
#> $error
#> [1] 9.976201e-05
#> 
#> $functionEvaluations
#> [1] 58377
#> 
#> $returnCode
#> [1] 0
#> 
pcubature(testFn5, rep(0,3), rep(1,3), tol=1e-4)
#> $integral
#> [1] 0.9999963
#> 
#> $error
#> [1] 7.175992e-05
#> 
#> $functionEvaluations
#> [1] 35937
#> 
#> $returnCode
#> [1] 0
#> 

##
## Test function 6 (Tsuda's example)
## Compare with original cubature result of
## ./cubature_test 4 1e-4 6 0
## 4-dim integral, tolerance = 0.0001
## integrand 6: integral = 0.999998, est err = 9.99685e-05, true err = 1.5717e-06
## #evals = 18753

testFn6 <- function(x) {
  a <- (1 + sqrt(10.0)) / 9.0
  prod(a / (a + 1) * ((a + 1) / (a + x))^2)
}

hcubature(testFn6, rep(0,4), rep(1,4), tol=1e-4)
#> $integral
#> [1] 0.9999984
#> 
#> $error
#> [1] 9.996851e-05
#> 
#> $functionEvaluations
#> [1] 18753
#> 
#> $returnCode
#> [1] 0
#> 
pcubature(testFn6, rep(0,4), rep(1,4), tol=1e-4)
#> Warning: pcubature not recommended for dimensions > 3!
#> $integral
#> [1] 1
#> 
#> $error
#> [1] 1.729307e-06
#> 
#> $functionEvaluations
#> [1] 83521
#> 
#> $returnCode
#> [1] 0
#> 


##
## Test function 7
##   test integrand from W. J. Morokoff and R. E. Caflisch, "Quasi=
##   Monte Carlo integration," J. Comput. Phys 122, 218-230 (1995).
##   Designed for integration on [0,1]^dim, integral = 1. */
## Compare with original cubature result of
## ./cubature_test 3 1e-4 7 0
## 3-dim integral, tolerance = 0.0001
## integrand 7: integral = 1.00001, est err = 9.96657e-05, true err = 1.15994e-05
## #evals = 7887

testFn7 <- function(x) {
  n <- length(x)
  p <- 1/n
  (1 + p)^n * prod(x^p)
}

hcubature(testFn7, rep(0,3), rep(1,3), tol=1e-4)
#> $integral
#> [1] 1.000012
#> 
#> $error
#> [1] 9.966567e-05
#> 
#> $functionEvaluations
#> [1] 7887
#> 
#> $returnCode
#> [1] 0
#> 
pcubature(testFn7, rep(0,3), rep(1,3), tol=1e-4)
#> $integral
#> [1] 0.9999878
#> 
#> $error
#> [1] 2.175551e-05
#> 
#> $functionEvaluations
#> [1] 274625
#> 
#> $returnCode
#> [1] 0
#> 


## Example from web page
## https://github.com/stevengj/cubature
##
## f(x) = exp(-0.5(euclidean_norm(x)^2)) over the three-dimensional
## hyperbcube [-2, 2]^3
## Compare with original cubature result
testFnWeb <-  function(x) {
  exp(-0.5 * sum(x^2))
}

hcubature(testFnWeb, rep(-2,3), rep(2,3), tol=1e-4)
#> $integral
#> [1] 13.69609
#> 
#> $error
#> [1] 0.001369187
#> 
#> $functionEvaluations
#> [1] 3399
#> 
#> $returnCode
#> [1] 0
#> 
pcubature(testFnWeb, rep(-2,3), rep(2,3), tol=1e-4)
#> $integral
#> [1] 13.69611
#> 
#> $error
#> [1] 7.242129e-05
#> 
#> $functionEvaluations
#> [1] 4913
#> 
#> $returnCode
#> [1] 0
#> 

## Test function I.1d from
## Numerical integration using Wang-Landau sampling
## Y. W. Li, T. Wust, D. P. Landau, H. Q. Lin
## Computer Physics Communications, 2007, 524-529
## Compare with exact answer: 1.63564436296
##
I.1d <- function(x) {
  sin(4*x) *
    x * ((x * ( x * (x*x-4) + 1) - 1))
}

hcubature(I.1d, -2, 2, tol=1e-7)
#> $integral
#> [1] 1.635644
#> 
#> $error
#> [1] 4.024021e-09
#> 
#> $functionEvaluations
#> [1] 105
#> 
#> $returnCode
#> [1] 0
#> 
pcubature(I.1d, -2, 2, tol=1e-7)
#> $integral
#> [1] 1.635644
#> 
#> $error
#> [1] 1.332268e-15
#> 
#> $functionEvaluations
#> [1] 65
#> 
#> $returnCode
#> [1] 0
#> 

## Test function I.2d from
## Numerical integration using Wang-Landau sampling
## Y. W. Li, T. Wust, D. P. Landau, H. Q. Lin
## Computer Physics Communications, 2007, 524-529
## Compare with exact answer: -0.01797992646
##
##
I.2d <- function(x) {
  x1 = x[1]
  x2 = x[2]
  sin(4*x1+1) * cos(4*x2) * x1 * (x1*(x1*x1)^2 - x2*(x2*x2 - x1) +2)
}

hcubature(I.2d, rep(-1, 2), rep(1, 2), maxEval=10000)
#> Warning: hcubature: maxEval (10000) exhausted before the requested tolerance was reached. Returning best-effort estimate; treat the result with caution and consider a larger maxEval, robust = TRUE, or cubintegrate(method = "cuhre").
#> $integral
#> [1] -0.01797993
#> 
#> $error
#> [1] 7.845607e-07
#> 
#> $functionEvaluations
#> [1] 10013
#> 
#> $returnCode
#> [1] 2
#> 
pcubature(I.2d, rep(-1, 2), rep(1, 2), maxEval=10000)
#> $integral
#> [1] -0.01797993
#> 
#> $error
#> [1] 2.092601e-10
#> 
#> $functionEvaluations
#> [1] 1089
#> 
#> $returnCode
#> [1] 0
#> 

##
## Example of multivariate normal integration borrowed from
## package mvtnorm (on CRAN) to check ... argument
## Compare with output of
## pmvnorm(lower=rep(-0.5, m), upper=c(1,4,2), mean=rep(0, m), corr=sigma, alg=Miwa())
##     0.3341125.  Blazing quick as well!  Ours is, not unexpectedly, much slower.
##
dmvnorm <- function (x, mean, sigma, log = FALSE) {
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    distval <- mahalanobis(x, center = mean, cov = sigma)
    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    if (log)
        return(logretval)
    exp(logretval)
}

m <- 3
sigma <- diag(3)
sigma[2,1] <- sigma[1, 2] <- 3/5 ; sigma[3,1] <- sigma[1, 3] <- 1/3
sigma[3,2] <- sigma[2, 3] <- 11/15
hcubature(dmvnorm, lower=rep(-0.5, m), upper=c(1,4,2),
                        mean=rep(0, m), sigma=sigma, log=FALSE,
               maxEval=10000)
#> Warning: hcubature: maxEval (10000) exhausted before the requested tolerance was reached. Returning best-effort estimate; treat the result with caution and consider a larger maxEval, robust = TRUE, or cubintegrate(method = "cuhre").
#> $integral
#> [1] 0.3341125
#> 
#> $error
#> [1] 4.185435e-06
#> 
#> $functionEvaluations
#> [1] 10065
#> 
#> $returnCode
#> [1] 2
#> 
pcubature(dmvnorm, lower=rep(-0.5, m), upper=c(1,4,2),
                        mean=rep(0, m), sigma=sigma, log=FALSE,
               maxEval=10000)
#> $integral
#> [1] 0.3341125
#> 
#> $error
#> [1] 2.21207e-06
#> 
#> $functionEvaluations
#> [1] 4913
#> 
#> $returnCode
#> [1] 0
#> 
# }
```
