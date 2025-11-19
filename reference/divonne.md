# Integration by Stratified Sampling for Variance Reduction

Divonne works by stratified sampling, where the partioning of the
integration region is aided by methods from numerical optimization.

## Usage

``` r
divonne(
  f,
  nComp = 1L,
  lowerLimit,
  upperLimit,
  ...,
  relTol = 1e-05,
  absTol = 1e-12,
  minEval = 0L,
  maxEval = 10^6,
  flags = list(verbose = 0L, final = 1L, keep_state = 0L, level = 0L),
  rngSeed = 0L,
  nVec = 1L,
  key1 = 47L,
  key2 = 1L,
  key3 = 1L,
  maxPass = 5L,
  border = 0,
  maxChisq = 10,
  minDeviation = 0.25,
  xGiven = NULL,
  nExtra = 0L,
  peakFinder = NULL,
  stateFile = NULL
)
```

## Arguments

- f:

  The function (integrand) to be integrated as in [`cuhre()`](cuhre.md).
  Optionally, the function can take an additional argument in addition
  to the variable being integrated: - `cuba_phase` - indicating the
  integration phase:

  0

  :   sampling of the points in `xgiven`

  1

  :   partitioning phase

  2

  :   final integration phase

  3

  :   refinement phase

  This information might be useful if the integrand takes long to
  compute and a sufficiently accurate approximation of the integrand is
  available. The actual value of the integral is only of minor
  importance in the partitioning phase, which is instead much more
  dependent on the peak structure of the integrand to find an
  appropriate tessellation. An approximation which reproduces the peak
  structure while leaving out the fine details might hence be a
  perfectly viable and much faster substitute when `cuba_phase < 2`. In
  all other instances, phase can be ignored and it is entirely
  admissible to define the integrand without it.

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

- key1:

  integer that determines sampling in the partitioning phase:
  `key1 = 7, 9, 11, 13` selects the cubature rule of degree `key1`. Note
  that the degree-11 rule is available only in 3 dimensions, the
  degree-13 rule only in 2 dimensions. For other values of `key1`, a
  quasi-random sample of \\n=\|key1\|\\ points is used, where the sign
  of `key1` determines the type of sample, `key1 = 0`, use the default
  rule. `key1 > 0`, use a Korobov quasi-random sample, `key1 < 0`, use a
  Sobol quasi-random sample if `flags$seed` is zero, otherwise a
  “standard” sample (Mersenne Twister) pseudo-random sample

- key2:

  integer that determines sampling in the final integration phase: same
  as `key1`, but here \\n=\|key2\|\\ determines the number of points,
  \\n \> 39\\, sample \\n\\ points, \\n \< 40\\, sample \\n\\ `nneed`
  points, where `nneed` is the number of points needed to reach the
  prescribed accuracy, as estimated by Divonne from the results of the
  partitioning phase.

- key3:

  integer that sets the strategy for the refinement phase: `key3 = 0`,
  do not treat the subregion any further. `key3 = 1`, split the
  subregion up once more. Otherwise, the subregion is sampled a third
  time with `key3` specifying the sampling parameters exactly as `key2`
  above.

- maxPass:

  integer that controls the thoroughness of the partitioning phase: The
  partitioning phase terminates when the estimated total number of
  integrand evaluations (partitioning plus final integration) does not
  decrease for `maxPass` successive iterations. A decrease in points
  generally indicates that Divonne discovered new structures of the
  integrand and was able to find a more effective partitioning.
  `maxPass` can be understood as the number of “safety” iterations that
  are performed before the partition is accepted as final and counting
  consequently restarts at zero whenever new structures are found.

- border:

  the relative width of the border of the integration region. Points
  falling into the border region will not be sampled directly, but will
  be extrapolated from two samples from the interior. Use a non-zero
  `border` if the integrand subroutine cannot produce values directly on
  the integration boundary. The relative width of the border is
  identical in all the dimensions. For example, set `border=0.1` for a
  border of width equal to 10\\ width of the integration region.

- maxChisq:

  the maximum \\\chi^2\\ value a single subregion is allowed to have in
  the final integration phase. Regions which fail this \\\chi^2\\ test
  and whose sample averages differ by more than `min.deviation` move on
  to the refinement phase.

- minDeviation:

  a bound, given as the fraction of the requested error of the entire
  integral, which determines whether it is worthwhile further examining
  a region that failed the \\\chi^2\\ test. Only if the two sampling
  averages obtained for the region differ by more than this bound is the
  region further treated.

- xGiven:

  a matrix (`nDim`, `nGiven`). A list of `nGiven` points where the
  integrand might have peaks. Divonne will consider these points when
  partitioning the integration region. The idea here is to help the
  integrator find the extrema of the integrand in the presence of very
  narrow peaks. Even if only the approximate location of such peaks is
  known, this can considerably speed up convergence.

- nExtra:

  the maximum number of extra points the peak-finder subroutine will
  return. If `nextra` is zero, `peakfinder` is not called and an
  arbitrary object may be passed in its place, e.g. just 0.

- peakFinder:

  the peak-finder subroutine. This R function is called whenever a
  region is up for subdivision and is supposed to point out possible
  peaks lying in the region, thus acting as the dynamic counterpart of
  the static list of points supplied in `xgiven`. It is expected to be
  declared as `peakFinder <- function(bounds, nMax)` where `bounds` is a
  matrix of dimension `(2, nDim)` which contains the lower (row 1) and
  upper (row 2) bounds of the subregion. The returned value should be a
  matrix `(nX, nDim)` where `nX` is the actual number of points (should
  be less or equal to `nMax`).

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

Divonne uses stratified sampling for variance reduction, that is, it
partitions the integration region such that all subregions have an
approximately equal value of a quantity called the spread (volume times
half-range).

See details in the documentation.

## References

J. H. Friedman, M. H. Wright (1981) A nested partitioning procedure for
numerical multiple integration. *ACM Trans. Math. Software*, **7**(1),
76-92.

J. H. Friedman, M. H. Wright (1981) User's guide for DIVONNE. SLAC
Report CGTM-193-REV, CGTM-193, Stanford University.

T. Hahn (2005) CUBA-a library for multidimensional numerical
integration. *Computer Physics Communications*, **168**, 78-95.

## See also

[`cuhre()`](cuhre.md), [`suave()`](suave.md), [`vegas()`](vegas.md)

## Examples

``` r
integrand <- function(arg, phase) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
}
divonne(integrand, relTol=1e-3,  absTol=1e-12, lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
        flags=list(verbose = 2),  key1= 47)
#> Divonne input parameters:
#>   ndim 3
#>   ncomp 1
#>   nvec 1
#>   epsrel 0.001
#>   epsabs 1e-12
#>   flags 6
#>   seed 0
#>   mineval 0
#>   maxeval 1000000
#>   key1 47
#>   key2 1
#>   key3 1
#>   maxpass 5
#>   border 0
#>   maxchisq 10
#>   mindeviation 0.25
#>   ngiven 0
#>   nextra 0
#>   statefile "(null)"
#> 
#> Partitioning phase:
#> 
#> Iteration 1 (pass 0):  8 regions
#>     836 integrand evaluations so far,
#>     406 in optimizing regions,
#>      70 in finding cuts
#> [1] 0.665011 +- 0.00470198
#> 
#> Iteration 2 (pass 0):  9 regions
#>     966 integrand evaluations so far,
#>     478 in optimizing regions,
#>      80 in finding cuts
#> [1] 0.664964 +- 0.00429467
#> 
#> Iteration 3 (pass 1):  10 regions
#>    1096 integrand evaluations so far,
#>     550 in optimizing regions,
#>      90 in finding cuts
#> [1] 0.664949 +- 0.00388393
#> 
#> Iteration 4 (pass 2):  11 regions
#>    1194 integrand evaluations so far,
#>     590 in optimizing regions,
#>     100 in finding cuts
#> [1] 0.664887 +- 0.0035842
#> 
#> Iteration 5 (pass 3):  12 regions
#>    1308 integrand evaluations so far,
#>     646 in optimizing regions,
#>     110 in finding cuts
#> [1] 0.664853 +- 0.0033385
#> 
#> Iteration 6 (pass 4):  13 regions
#>    1438 integrand evaluations so far,
#>     718 in optimizing regions,
#>     120 in finding cuts
#> [1] 0.664825 +- 0.0031276
#> 
#> Iteration 7 (pass 5):  14 regions
#>    1568 integrand evaluations so far,
#>     790 in optimizing regions,
#>     130 in finding cuts
#> [1] 0.664816 +- 0.00292074
#> 
#> Main integration on 14 regions with 211 samples per region.
#> $integral
#> [1] 0.6646098
#> 
#> $error
#> [1] 0.0006505917
#> 
#> $neval
#> [1] 3052
#> 
#> $prob
#> [1] 1.110223e-16
#> 
#> $returnCode
#> [1] 0
#> 

# Example with a peak-finder function
nDim <- 3L
peakf <- function(bounds, nMax) {
#  print(bounds) # matrix (ndim,2)
  x <- matrix(0, ncol = nMax, nrow = nDim)
   pas <- 1 / (nMax - 1)
   # 1ier point
   x[, 1] <- rep(0, nDim)
   # Les autres points
   for (i in 2L:nMax) {
      x[, i] <- x[, (i - 1)] + pas
    }
  x
} #end peakf

divonne(integrand, relTol=1e-3,  absTol=1e-12,
        lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
        flags=list(verbose = 2),  peakFinder = peakf, nExtra = 4L)
#> Divonne input parameters:
#>   ndim 3
#>   ncomp 1
#>   nvec 1
#>   epsrel 0.001
#>   epsabs 1e-12
#>   flags 6
#>   seed 0
#>   mineval 0
#>   maxeval 1000000
#>   key1 47
#>   key2 1
#>   key3 1
#>   maxpass 5
#>   border 0
#>   maxchisq 10
#>   mindeviation 0.25
#>   ngiven 0
#>   nextra 4
#>   statefile "(null)"
#> 
#> Partitioning phase:
#> 
#> Iteration 1 (pass 0):  8 regions
#>     857 integrand evaluations so far,
#>     382 in optimizing regions,
#>      70 in finding cuts
#> [1] 0.665011 +- 0.00470198
#> 
#> Iteration 2 (pass 0):  9 regions
#>     993 integrand evaluations so far,
#>     454 in optimizing regions,
#>      80 in finding cuts
#> [1] 0.664964 +- 0.00429467
#> 
#> Iteration 3 (pass 1):  10 regions
#>    1129 integrand evaluations so far,
#>     526 in optimizing regions,
#>      90 in finding cuts
#> [1] 0.664949 +- 0.00388393
#> 
#> Iteration 4 (pass 2):  11 regions
#>    1233 integrand evaluations so far,
#>     566 in optimizing regions,
#>     100 in finding cuts
#> [1] 0.664887 +- 0.0035842
#> 
#> Iteration 5 (pass 3):  12 regions
#>    1353 integrand evaluations so far,
#>     622 in optimizing regions,
#>     110 in finding cuts
#> [1] 0.664853 +- 0.0033385
#> 
#> Iteration 6 (pass 4):  13 regions
#>    1489 integrand evaluations so far,
#>     694 in optimizing regions,
#>     120 in finding cuts
#> [1] 0.664825 +- 0.0031276
#> 
#> Iteration 7 (pass 5):  14 regions
#>    1625 integrand evaluations so far,
#>     766 in optimizing regions,
#>     130 in finding cuts
#> [1] 0.664816 +- 0.00292074
#> 
#> Main integration on 14 regions with 211 samples per region.
#> $integral
#> [1] 0.6646098
#> 
#> $error
#> [1] 0.0006505917
#> 
#> $neval
#> [1] 3109
#> 
#> $prob
#> [1] 1.110223e-16
#> 
#> $returnCode
#> [1] 0
#> 
```
