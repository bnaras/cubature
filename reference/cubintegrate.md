# Unified Cubature Integration Interface

Integrate a function within specified limits using method specified.
Further arguments specific to method as well as other arguments to f may
be passed. For defaults used in each method, see help on the method or
[`default_args()`](default_args.md).

## Usage

``` r
cubintegrate(
  f,
  lower,
  upper,
  fDim = 1,
  method = c("hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas"),
  relTol = 1e-05,
  absTol = 1e-12,
  maxEval = 10^6,
  nVec = 1L,
  ...
)
```

## Arguments

- f:

  The function (integrand) to be integrated. Can be vectorized version,
  but the additional arguments `...` must indicate via either
  `vectorInterface = TRUE` for `hcubature` and `pcubature`, or a value
  for `nVec`. See details on each method.

- lower:

  The lower limit of integration, a vector for hypercubes.

- upper:

  The upper limit of integration, a vector for hypercubes.

- fDim:

  The number of components of f, default 1, bears no relation to the
  dimension of the hypercube over which integration is performed.

- method:

  the method to use should be one of "hcubature", "pcubature", "cuhre",
  "divonne", "suave" or "vegas".

- relTol:

  The maximum tolerance, default 1e-5.

- absTol:

  the absolute tolerance, default 1e-12.

- maxEval:

  The maximum number of function evaluations needed, default 10^6. Note
  that the actual number of function evaluations performed is only
  approximately guaranteed not to exceed this number.

- nVec:

  the number of vectorization points for Cuba C library, default 1, but
  can be set to an integer \> 1 for vectorization, for example, 1024.
  The function f above needs to handle the vector of points
  appropriately; see vignette examples. Unlike Cuba, the cubature C
  library manages the number of points on its own and can vary between
  calls. Therefore, any value for nVec greater than one implies
  vectorization for a cubature method.

- ...:

  All other arguments which may include integration method specific
  parameters and those for f. Unrecognized parameters for integration
  method are presumed to be intended for f and so processed.

## Value

The returned value is a list of items:

- integral:

  the value of the integral

- error:

  the estimated absolute error

- neval:

  the number of times the function was evaluated

- returnCode:

  the actual integer return code of the C routine; a non-zero value
  usually indicates problems; further interpretation depends on method

- nregions:

  forcCuba routines, the actual number of subregions needed

- prob:

  the \\\chi^2\\-probability (not the \\\chi^2\\-value itself!) that
  `error` is not a reliable estimate of the true integration error.

## See also

[`default_args()`](default_args.md), [`hcubature()`](hcubature.md),
[`pcubature()`](hcubature.md), [`cuhre()`](cuhre.md),
[`vegas()`](vegas.md), [`suave()`](suave.md), [`divonne()`](divonne.md)

## Examples

``` r
I.1d <- function(x) {
  sin(4*x) *
    x * ((x * ( x * (x*x-4) + 1) - 1))
}
I.1d_v <- function(x) {
   matrix(apply(x, 2, function(z)
       sin(4 * z) *
       z * ((z * ( z * (z * z - 4) + 1) - 1))),
       ncol = ncol(x))
}
cubintegrate(f = I.1d, lower = -2, upper = 2, method = "pcubature")
#> $integral
#> [1] 1.635644
#> 
#> $error
#> [1] 1.332268e-15
#> 
#> $neval
#> [1] 65
#> 
#> $returnCode
#> [1] 0
#> 
cubintegrate(f = I.1d, lower = -2, upper = 2, method = "cuhre", flags=list(verbose = 2))
#> Cuhre input parameters:
#>   ndim 1
#>   ncomp 1
#>   nvec 1
#>   epsrel 1e-05
#>   epsabs 1e-12
#>   flags 6
#>   mineval 0
#>   maxeval 1000000
#>   key 0
#>   statefile "(null)"
#> 
#> Iteration 1:  11 integrand evaluations so far
#> [1] 1.79496 +- 36.4275   chisq 0 (0 df)
#> 
#> Iteration 2:  33 integrand evaluations so far
#> [1] 1.60939 +- 15.0794   chisq 2.21553e-05 (1 df)
#> 
#> Iteration 3:  55 integrand evaluations so far
#> [1] 1.62346 +- 7.01714   chisq 2.28586e-05 (2 df)
#> 
#> Iteration 4:  77 integrand evaluations so far
#> [1] 1.63569 +- 0.145134      chisq 2.5197e-05 (3 df)
#> 
#> Iteration 5:  99 integrand evaluations so far
#> [1] 1.63568 +- 0.0643672     chisq 2.51973e-05 (4 df)
#> 
#> Iteration 6:  121 integrand evaluations so far
#> [1] 1.63568 +- 0.0123924     chisq 2.52128e-05 (5 df)
#> 
#> Iteration 7:  143 integrand evaluations so far
#> [1] 1.63566 +- 0.00587027    chisq 2.67284e-05 (6 df)
#> 
#> Iteration 8:  165 integrand evaluations so far
#> [1] 1.63564 +- 5.32936e-05   chisq 3.85791e-05 (7 df)
#> 
#> Iteration 9:  187 integrand evaluations so far
#> [1] 1.63564 +- 3.85904e-05   chisq 3.8594e-05 (8 df)
#> 
#> Iteration 10:  209 integrand evaluations so far
#> [1] 1.63564 +- 2.5705e-05    chisq 3.86238e-05 (9 df)
#> 
#> Iteration 11:  231 integrand evaluations so far
#> [1] 1.63564 +- 2.02465e-05   chisq 3.86238e-05 (10 df)
#> 
#> Iteration 12:  253 integrand evaluations so far
#> [1] 1.63564 +- 1.51801e-05   chisq 3.8147e-05 (11 df)
#> $integral
#> [1] 1.635644
#> 
#> $error
#> [1] 1.518009e-05
#> 
#> $nregions
#> [1] 12
#> 
#> $neval
#> [1] 253
#> 
#> $prob
#> [1] 0
#> 
#> $returnCode
#> [1] 0
#> 
cubintegrate(f = I.1d_v, lower = -2, upper = 2, method = "hcubature", nVec = 2L)
#> $integral
#> [1] 1.635644
#> 
#> $error
#> [1] 4.024021e-09
#> 
#> $neval
#> [1] 105
#> 
#> $returnCode
#> [1] 0
#> 
cubintegrate(f = I.1d_v, lower = -2, upper = 2, method = "cuhre", nVec = 128L)
#> $integral
#> [1] 1.635644
#> 
#> $error
#> [1] 1.518009e-05
#> 
#> $nregions
#> [1] 12
#> 
#> $neval
#> [1] 253
#> 
#> $prob
#> [1] 0
#> 
#> $returnCode
#> [1] 0
#> 
```
