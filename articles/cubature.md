# Get started with cubature

## Introduction

`cubature` approximates definite integrals over hypercubes in arbitrary
dimension, including improper integrals over infinite limits, by
evaluating the integrand at adaptively-chosen sample points. It bundles
two well-known C libraries:

- the [cubature library](https://github.com/stevengj/cubature) of
  Steven G. Johnson (routines `hcubature` and `pcubature`), and
- the [Cuba library](https://feynarts.de/cuba/) of Thomas Hahn (routines
  `cuhre`, `divonne`, `suave`, `vegas`).

Together these give you six integration methods in a single package.
Each method is available through its own named function (`hcubature`,
`pcubature`, `cuhre`, `divonne`, `suave`, `vegas`) as well as through a
uniform dispatcher, [`cubintegrate()`](../reference/cubintegrate.md),
that lets you switch methods with a single argument. Vectorized variants
of every method are also available for performance-sensitive work; see
the [vectorization vignette](vectorization.md) for timing data and
recommendations.

## Your first integral

The simplest possible call integrates a one-argument function over an
interval.

``` r
library(cubature)
hcubature(f = sin, lowerLimit = 0, upperLimit = pi)
```

    ## $integral
    ## [1] 2
    ## 
    ## $error
    ## [1] 2.220446e-14
    ## 
    ## $functionEvaluations
    ## [1] 15
    ## 
    ## $returnCode
    ## [1] 0

The returned list has four elements:

- **`integral`** — the estimated value of the integral,
- **`error`** — the estimated absolute error,
- **`functionEvaluations`** — how many times the integrand was called,
- **`returnCode`** — an integer status code, where `0` means the
  algorithm reached the requested tolerance. A non-zero value means
  something went wrong; see [`?hcubature`](../reference/hcubature.md)
  for the full meaning (`1` = hard internal failure, `2` = `maxEval`
  budget exhausted before tolerance was met).

Here is a two-dimensional example that integrates the product
`cos(x₁)·cos(x₂)` over the unit square.

``` r
f <- function(x) cos(x[1]) * cos(x[2])
hcubature(f, lowerLimit = c(0, 0), upperLimit = c(1, 1))
```

    ## $integral
    ## [1] 0.7080734
    ## 
    ## $error
    ## [1] 2.888228e-06
    ## 
    ## $functionEvaluations
    ## [1] 51
    ## 
    ## $returnCode
    ## [1] 0

The exact answer is $\sin(1)^{2} \approx 0.7081$.

The integration limits can be infinite — `hcubature` and `pcubature`
automatically transform the variable using `atan` so the unit-hypercube
routines inside still work:

``` r
## Standard normal density integrates to 1 over the real line
hcubature(function(x) exp(-x^2 / 2) / sqrt(2 * pi),
          lowerLimit = -Inf, upperLimit = Inf)
```

    ## $integral
    ## [1] 1
    ## 
    ## $error
    ## [1] 7.709821e-08
    ## 
    ## $functionEvaluations
    ## [1] 225
    ## 
    ## $returnCode
    ## [1] 0

## The unified `cubintegrate()` interface

The package exposes a second, uniform way to call every integration
method: [`cubintegrate()`](../reference/cubintegrate.md). One signature,
one argument convention, and a `method =` argument to pick which C
routine you want. This is the interface we recommend for new code,
because swapping methods to compare results or tune performance becomes
a one-argument change.

``` r
## The same 2D integral via three different methods, all the same way
cubintegrate(f, lower = c(0, 0), upper = c(1, 1), method = "hcubature")$integral
```

    ## [1] 0.7080734

``` r
cubintegrate(f, lower = c(0, 0), upper = c(1, 1), method = "cuhre")$integral
```

    ## [1] 0.7080734

``` r
cubintegrate(f, lower = c(0, 0), upper = c(1, 1), method = "pcubature")$integral
```

    ## [1] 0.7080734

Each of those would otherwise have required `hcubature(...)`,
`cuhre(...)`, `pcubature(...)` — different function names, slightly
different argument names, and (for the Cuba routines) different
conventions for passing method-specific options. `cubintegrate` papers
over all of that.

Method-specific arguments (e.g. Cuhre’s `key` parameter or Divonne’s
`maxPass`) are passed through `...`. To see what each method accepts,
use [`default_args()`](../reference/default_args.md):

``` r
str(default_args(), list.len = 4, vec.len = 2)
```

    ## List of 6
    ##  $ hcubature:List of 1
    ##   ..$ norm: chr [1:5] "INDIVIDUAL" "PAIRED" ...
    ##  $ pcubature:List of 1
    ##   ..$ norm: chr [1:5] "INDIVIDUAL" "PAIRED" ...
    ##  $ cuhre    :List of 4
    ##   ..$ minEval  : int 0
    ##   ..$ stateFile: NULL
    ##   ..$ flags    :List of 6
    ##   .. ..$ verbose   : int 0
    ##   .. ..$ final     : int 1
    ##   .. ..$ smooth    : int 0
    ##   .. ..$ keep_state: int 0
    ##   .. .. [list output truncated]
    ##   ..$ key      : int 0
    ##  $ divonne  :List of 14
    ##   ..$ minEval     : int 0
    ##   ..$ stateFile   : NULL
    ##   ..$ flags       :List of 6
    ##   .. ..$ verbose   : int 0
    ##   .. ..$ final     : int 1
    ##   .. ..$ smooth    : int 0
    ##   .. ..$ keep_state: int 0
    ##   .. .. [list output truncated]
    ##   ..$ rngSeed     : int 0
    ##   .. [list output truncated]
    ##   [list output truncated]

## Choosing a method

The six methods are not interchangeable — each is best on a different
class of integrand. A short decision table:

| Situation                             | Recommended                                 |
|---------------------------------------|---------------------------------------------|
| You’re not sure / general use         | `hcubature`                                 |
| Smooth integrand, low dimension (≤ 3) | `pcubature` or `cuhre`                      |
| Non-smooth, localized, or suspicious  | `cuhre`, or `hcubature(..., robust = TRUE)` |
| High dimension (\> 5)                 | `vegas`, `suave`, `divonne` (Monte Carlo)   |
| You want a second opinion             | any of the above, compared                  |

`hcubature` uses h-adaptive Genz–Malik subdivision; it refines the
integration region where the integrand varies most. `pcubature`
increases the polynomial degree of a tensor Clenshaw–Curtis rule on a
fixed domain; it can be very efficient on smooth integrands but is
brittle when features are localized (see the robustness section below).
Cuba’s `cuhre` is h-adaptive like `hcubature` but with a denser base
rule, which makes it substantially more robust on non-smooth integrands
— at the cost of more function evaluations per subregion. The Monte
Carlo methods (`vegas`, `suave`, `divonne`) scale better than
deterministic methods in high dimensions but return probabilistic error
estimates rather than deterministic bounds.

When in genuine doubt, run your integrand through `cubintegrate` with
two different methods and check they agree.

## Vector-valued integrands

All six methods support integrands that return more than one scalar per
call — `fDim > 1` for `hcubature`/`pcubature`, `nComp > 1` for the Cuba
methods. The integrand should return a numeric vector of length `fDim`
(scalar interface) or a `fDim × npt` matrix (vectorized interface).
Example:

``` r
## Simultaneously integrate sin and cos over [0, pi/2]
f2 <- function(x) c(sin(x), cos(x))
hcubature(f2, 0, pi / 2, fDim = 2)$integral
```

    ## [1] 1 1

``` r
## Exact: both equal 1
```

## Robustness: when hcubature and pcubature can silently fail

The default error estimator in `hcubature` uses a single “null rule”
difference between the degree-7 Genz–Malik rule and its embedded
degree-5 companion on each region. For smooth integrands this is
reliable and cheap. For integrands with **localized support** — a thin
hinge, a near-step sigmoid, a very narrow Gaussian inside an oversized
bounding box — it can fail in a particular way: the rule’s sample points
miss the integrand’s support entirely, all sample values are zero, the
rule estimate is zero with zero estimated error, and the algorithm
confidently reports “converged” on a completely wrong answer. This is
distinct from ordinary numerical error — the algorithm is internally
consistent; it has just never evaluated the function at a point where
it’s nonzero. `pcubature` can fail the same way because its
Clenshaw–Curtis tensor grid is also a fixed sample geometry.

A small reproducible example: a hinge-squared times a steep logistic
times a standard Gaussian over `[-10, 10]²`. With default settings
`hcubature` returns exactly `0` (with `returnCode = 0`, i.e. silent
“success”), while the true integral is about 0.21:

``` r
f <- function(x) {
  t1 <- -3.5862; t2 <- -9.3801
  beta0 <- 9.5;  gamma0 <- 19 * sqrt(0.75)
  (max(t1 + t2 * x[1] - x[2], 0))^2 *
    plogis(beta0 + gamma0 * x[1]) *
    dnorm(x[1]) * dnorm(x[2])
}
## Default: silently wrong
hcubature(f, c(-10, -10), c(10, 10))$integral
#> [1] 0
## robust = TRUE: correct
hcubature(f, c(-10, -10), c(10, 10), robust = TRUE)$integral
#> [1] ~0.21
## cuhre: also correct
cubintegrate(f, c(-10, -10), c(10, 10), method = "cuhre")$integral
#> [1] ~0.21
```

### The opt-in `robust` path

Starting in cubature 2.2.0, [`hcubature()`](../reference/hcubature.md),
[`pcubature()`](../reference/hcubature.md), and
[`cubintegrate()`](../reference/cubintegrate.md) accept an optional
`robust = TRUE` argument that activates extra safeguards. For
`hcubature` these are: an additional degree-3 null rule combined with
the existing degree-5 one via a Cuhre-style decay-check formula
(Berntsen, Espelid & Genz, *ACM TOMS* **17**(4), 437–451, 1991); a
parent–child consistency check at each subdivision step that inflates
children’s error estimates by the discrepancy between their sum and the
parent’s prior estimate; and a running peak-density “suspect region”
detector. For `pcubature` the robust path raises the initial
Clenshaw–Curtis tensor order so the starting grid is denser. In both
cases the default (non-robust) path is unchanged bit-for-bit; the
safeguards are paid for only when opted into.

The cost of `robust = TRUE` on smooth integrands is near zero in most
cases and at most ~25% extra function evaluations on sharply-peaked
ones. The cost on genuinely problem integrands (where the default path
was silently wrong) is several times more function evaluations — which
is the price of correctness. When `method` is a Cuba method,
`cubintegrate` silently ignores `robust` because those integrators
already use more conservative error estimators.

### When robust is not enough: the sampling-density cliff

Neither robust path is a complete cure. The safeguards work on regions
where the rule has *seen some signal* and needs a better estimate of it.
If the integrand’s support is so narrow that every sample point of the
base rule misses it entirely, the safeguards have nothing to work with.
This failure mode is determined by the rule’s sample geometry. For
`hcubature`’s Genz–Malik rule on `[-h, h]^d`, the cliff for hinge-line
integrands sits at `t₁ = -√(9/70)·h ≈ -0.3586·h`; for `pcubature`’s
robust (m = 2) Clenshaw–Curtis grid it’s at
`t₁ = -cos(3π/8)·h ≈ -0.3827·h`. For parameters on one side of the cliff
every sample point fires and the algorithm converges correctly; for
parameters on the other side every sample point evaluates to zero and
the algorithm silently reports zero. No combination of error-estimator
tweaks can detect this.

The only way to address the cliff is a **denser base rule**, which is
exactly what Cuba’s Cuhre provides: it uses a degree-13 rule with 65
sample points in 2D, far more than Genz–Malik’s 17, making the cliff
arbitrarily thin. **If your integrand is non-smooth and has features
much smaller than the integration box, prefer
`cubintegrate(method = "cuhre")` over `hcubature` / `pcubature` — even
with `robust = TRUE`.** As a secondary workaround, shrinking the
integration box (removing tail regions that contribute negligibly)
reduces the scale mismatch between the rule and the integrand’s features
and can make either `hcubature` or `pcubature` reliable on problems
where they would otherwise fail.

### The `NOT_CONVERGED` return code

Independently of `robust`, starting in 2.2.0 `hcubature` and `pcubature`
return `returnCode = 2` (and emit an R
[`warning()`](https://rdrr.io/r/base/warning.html)) when the `maxEval`
budget is exhausted without meeting the requested tolerance. Previously
such cases silently returned `returnCode = 0` (success), with the caller
expected to compare `error` against their own tolerance. The new
behavior is a strict improvement in honesty: the `error` field is still
populated with the best available estimate, and the returned `integral`
is still the best-effort result, so users who want the permissive old
behavior can wrap the call in `suppressWarnings(...)` and/or check
`returnCode %in% c(0, 2)`. But callers that defensively did
`if (returnCode == 0) use(result)` will now correctly fall through on
non-convergence rather than silently accept a possibly-wrong answer. The
symbolic constants `CUBATURE_SUCCESS`, `CUBATURE_FAILURE`, and
`CUBATURE_NOT_CONVERGED` are exported from the upstream C header for
compiled clients.

## Where to go next

- The [vectorization vignette](vectorization.md) benchmarks scalar vs.
  vectorized interfaces across ten test integrands and makes per-method
  recommendations — read this if you call the integrator many times in a
  loop and want maximum performance.
- The reference pages for
  \[[`hcubature()`](../reference/hcubature.md)\],
  \[[`pcubature()`](../reference/hcubature.md)\],
  \[[`cubintegrate()`](../reference/cubintegrate.md)\],
  \[[`cuhre()`](../reference/cuhre.md)\],
  \[[`divonne()`](../reference/divonne.md)\],
  \[[`suave()`](../reference/suave.md)\], and
  \[[`vegas()`](../reference/vegas.md)\] document every method-specific
  option.
- The source for both bundled C libraries is in `src/libcubature/` and
  `src/Cuba/`.

## Implementation notes

The only real modification we have made to the underlying `cubature`
library beyond the 2.2.0 robustness safeguards is that we use `M = 16`
rather than the default `M = 19` suggested by the original author for
`pcubature`, to comply with CRAN package size limits. The changes made
to the `Cuba` library are managed in a [GitHub repo
branch](https://github.com/bnaras/Cuba/tree/R_pkg): each time a new
release is made, we update the main branch and keep Unix-platform
changes in a branch named `R_pkg` against the current main branch.
Windows customization lives in `Makevars.win`.

## Session info

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] cubature_2.2.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    ##  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
    ##  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
    ## [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
    ## [17] compiler_4.5.3    tools_4.5.3       ragg_1.5.2        bslib_0.10.0     
    ## [21] evaluate_1.0.5    Rcpp_1.1.1        yaml_2.3.12       jsonlite_2.0.0   
    ## [25] rlang_1.2.0       fs_2.0.1
