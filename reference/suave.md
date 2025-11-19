# Integration with SUbregion-Adaptive Vegas Algorithm

Suave uses [`vegas()`](vegas.md)-like importance sampling combined with
a globally adaptive subdivision strategy: Until the requested accuracy
is reached, the region with the largest error at the time is bisected in
the dimension in which the fluctuations of the integrand are reduced
most. The number of new samples in each half is prorated for the
fluctuation in that half.

## Usage

``` r
suave(
  f,
  nComp = 1L,
  lowerLimit,
  upperLimit,
  ...,
  relTol = 1e-05,
  absTol = 1e-12,
  minEval = 0L,
  maxEval = 10^6,
  flags = list(verbose = 0L, final = 1L, smooth = 0L, keep_state = 0L, level = 0L),
  rngSeed = 0L,
  nVec = 1L,
  nNew = 1000L,
  nMin = 50L,
  flatness = 50,
  stateFile = NULL
)
```

## Arguments

- f:

  The function (integrand) to be integrated as in [`cuhre()`](cuhre.md).
  Optionally, the function can take two additional arguments in addition
  to the variable being integrated: - `cuba_weight` which is the weight
  of the point being sampled, - `cuba_iter` the current iteration
  number. The function author may choose to use these in any appropriate
  way or ignore them altogether.

- nComp:

  The number of components of f, default 1, bears no relation to the
  dimension of the hypercube over which integration is performed.

- lowerLimit:

  The lower limit of integration, a vector for hypercubes.

- upperLimit:

  The upper limit of integration, a vector for hypercubes.

- ...:

  All other arguments passed to the function f.

- relTol:

  The maximum tolerance, default 1e-5.

- absTol:

  the absolute tolerance, default 1e-12.

- minEval:

  the minimum number of function evaluations required

- maxEval:

  The maximum number of function evaluations needed, default 10^6. Note
  that the actual number of function evaluations performed is only
  approximately guaranteed not to exceed this number.

- flags:

  flags governing the integration. The list here is exhaustive to keep
  the documentation and invocation uniform, but not all flags may be
  used for a particular method as noted below. List components:

  verbose

  :   encodes the verbosity level, from 0 (default) to 3. Level 0 does
      not print any output, level 1 prints reasonable information on the
      progress of the integration, level 2 also echoes the input
      parameters, and level 3 further prints the subregion results.

  final

  :   when 0, all sets of samples collected on a subregion during the
      various iterations or phases contribute to the final result. When
      1, only the last (largest) set of samples is used in the final
      result.

  smooth

  :   Applies to Suave and Vegas only. When 0, apply additional
      smoothing to the importance function, this moderately improves
      convergence for many integrands. When 1, use the importance
      function without smoothing, this should be chosen if the integrand
      has sharp edges.

  keep_state

  :   when nonzero, retain state file if argument `stateFile` is
      non-null, else delete `stateFile` if specified.

  load_state

  :   Applies to Vegas only. Reset the integrator state even if a state
      file is present, i.e. keep only the grid. Together with
      `keep_state` this allows a grid adapted by one integration to be
      used for another integrand.

  level

  :   applies only to Divonne, Suave and Vegas. When `0`, Mersenne
      Twister random numbers are used. When nonzero Ranlux random
      numbers are used, except when `rngSeed` is zero which forces use
      of Sobol quasi-random numbers. Ranlux implements Marsaglia and
      Zaman's 24-bit RCARRY algorithm with generation period \\p\\, i.e.
      for every 24 generated numbers used, another \\p-24\\ are skipped.
      The luxury level for the Ranlux generator may be encoded in
      `level` as follows:

      Level 1 (p = 48)

      :   gives very long period, passes the gap test but fails spectral
          test

      Level 2 (p = 97)

      :   passes all known tests, but theoretically still defective

      Level 3 (p = 223)

      :   any theoretically possible correlations have very small chance
          of being observed

      Level 4 (p = 389)

      :   highest possible luxury, all 24 bits chaotic

      Levels 5-23

      :   default to 3, values above 24 directly specify the period p.
          Note that Ranlux's original level 0, (mis)used for selecting
          Mersenne Twister in Cuba, is equivalent to `level = 24`

- rngSeed:

  seed, default 0, for the random number generator. Note the
  articulation with `level` settings for `flag`

- nVec:

  the number of vectorization points, default 1, but can be set to an
  integer \> 1 for vectorization, for example, 1024 and the function f
  above needs to handle the vector of points appropriately. See vignette
  examples.

- nNew:

  the number of new integrand evaluations in each subdivision.

- nMin:

  the minimum number of samples a former pass must contribute to a
  subregion to be considered in that region's compound integral value.
  Increasing nmin may reduce jumps in the \\\chi^2\\ value.

- flatness:

  the parameter p, or the type of norm used to compute the fluctuation
  of a sample. This determines how prominently "outliers," i.e.
  individual samples with a large fluctuation, figure in the total
  fluctuation, which in turn determines how a region is split up. As
  suggested by its name, flatness should be chosen large for "flat"
  integrands and small for "volatile" integrands with high peaks. Note
  that since flatness appears in the exponent, one should not use too
  large values (say, no more than a few hundred) lest terms be truncated
  internally to prevent overflow.

- stateFile:

  the name of an external file. Vegas can store its entire internal
  state (i.e. all the information to resume an interrupted integration)
  in an external file. The state file is updated after every iteration.
  If, on a subsequent invocation, Vegas finds a file of the specified
  name, it loads the internal state and continues from the point it left
  off. Needless to say, using an existing state file with a different
  integrand generally leads to wrong results. Once the integration
  finishes successfully, i.e. the prescribed accuracy is attained, the
  state file is removed. This feature is useful mainly to define
  ‘check-points’ in long-running integrations from which the calculation
  can be restarted.

## Value

A list with components:

- neval:

  the actual number of integrand evaluations needed

- returnCode:

  if zero, the desired accuracy was reached, if -1, dimension out of
  range, if 1, the accuracy goal was not met within the allowed maximum
  number of integrand evaluations.

- integral:

  vector of length `nComp`; the integral of `integrand` over the
  hypercube

- error:

  vector of length `nComp`; the presumed absolute error of `integral`

- prob:

  vector of length `nComp`; the \\\chi^2\\-probability (not the
  \\\chi^2\\-value itself!) that `error` is not a reliable estimate of
  the true integration error.

## Details

See details in the documentation.

## References

T. Hahn (2005) CUBA-a library for multidimensional numerical
integration. *Computer Physics Communications*, **168**, 78-95.

## See also

[`cuhre()`](cuhre.md), [`divonne()`](divonne.md), [`vegas()`](vegas.md)

## Examples

``` r
integrand <- function(arg) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # end integrand
suave(integrand, lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
             relTol=1e-3,  absTol=1e-12,
             flags=list(verbose=2, final=0))
#> Suave input parameters:
#>   ndim 3
#>   ncomp 1
#>   nvec 1
#>   epsrel 0.001
#>   epsabs 1e-12
#>   flags 2
#>   seed 0
#>   mineval 0
#>   maxeval 1000000
#>   nnew 1000
#>   nmin 50
#>   flatness 50
#>   statefile "(null)"
#> 
#> Iteration 1:  1000 integrand evaluations so far
#> [1] 0.664916 +- 0.0138647    chisq 0 (0 df)
#> 
#> Iteration 2:  2000 integrand evaluations so far
#> [1] 0.664522 +- 0.00338364   chisq 0.00368404 (2 df)
#> 
#> Iteration 3:  3000 integrand evaluations so far
#> [1] 0.665025 +- 0.00238104   chisq 0.0739355 (5 df)
#> 
#> Iteration 4:  4000 integrand evaluations so far
#> [1] 0.664653 +- 0.0011072    chisq 0.228312 (8 df)
#> 
#> Iteration 5:  5000 integrand evaluations so far
#> [1] 0.664638 +- 0.000945872      chisq 1.56745 (12 df)
#> 
#> Iteration 6:  6000 integrand evaluations so far
#> [1] 0.66473 +- 0.000749366   chisq 5.46194 (16 df)
#> 
#> Iteration 7:  7000 integrand evaluations so far
#> [1] 0.66462 +- 0.000580976   chisq 14.8927 (20 df)
#> $integral
#> [1] 0.6646196
#> 
#> $error
#> [1] 0.0005809762
#> 
#> $neval
#> [1] 7000
#> 
#> $prob
#> [1] 0.2174863
#> 
#> $returnCode
#> [1] 0
#> 
```
