# Vectorization benchmarks

## When to read this

This vignette assumes you’ve already read the [Get started](cubature.md)
vignette and you now have a specific performance question: **is it worth
vectorizing my R integrand?**

`cubature` exposes two calling conventions for the same underlying C
routine:

- **Scalar interface** (the default) — your R function receives one
  point at a time and returns one scalar (or length-`fDim` vector).
- **Vectorized interface** (opt-in, `vectorInterface = TRUE` for
  `hcubature`/`pcubature` or `nVec > 1` for Cuba methods) — your R
  function receives a matrix of points and returns a matrix of values,
  amortizing R’s per-call overhead over many integrand evaluations at
  once.

The vectorized interface almost always wins for non-trivial problems,
sometimes by a factor of 10–50×, because the dominant cost is R’s
function-call overhead rather than the arithmetic inside the integrand.
But the speedup depends on your specific integrand, the method you pick,
and the dimension. This vignette runs ten representative integrands
through every combination and reports timings, so you can make an
informed decision rather than guess.

## A timing harness

Our harness will provide timing results for `hcubature`, `pcubature`
(where appropriate), and Cuba `cuhre` calls — scalar and vectorized
variants of each.

``` r
library(bench)
library(cubature)

harness <- function(which = NULL,
                    f, fv, lowerLimit, upperLimit, tol = 1e-3, iterations = 20, ...) {

  fns <- c(hc = "Non-vectorized Hcubature",
           hc.v = "Vectorized Hcubature",
           pc = "Non-vectorized Pcubature",
           pc.v = "Vectorized Pcubature",
           cc = "Non-vectorized cubature::cuhre",
           cc_v = "Vectorized cubature::cuhre")
  cc <- function() cubature::cuhre(f = f,
                                   lowerLimit = lowerLimit, upperLimit = upperLimit,
                                   relTol = tol,
                                   ...)
  cc_v <- function() cubature::cuhre(f = fv,
                                     lowerLimit = lowerLimit, upperLimit = upperLimit,
                                     relTol = tol,
                                     nVec = 1024L,
                                     ...)
  hc <- function() cubature::hcubature(f = f,
                                       lowerLimit = lowerLimit,
                                       upperLimit = upperLimit,
                                       tol = tol,
                                       ...)
  hc.v <- function() cubature::hcubature(f = fv,
                                         lowerLimit = lowerLimit,
                                         upperLimit = upperLimit,
                                         tol = tol,
                                         vectorInterface = TRUE,
                                         ...)
  pc <- function() cubature::pcubature(f = f,
                                       lowerLimit = lowerLimit,
                                       upperLimit = upperLimit,
                                       tol = tol,
                                       ...)
  pc.v <- function() cubature::pcubature(f = fv,
                                         lowerLimit = lowerLimit,
                                         upperLimit = upperLimit,
                                         tol = tol,
                                         vectorInterface = TRUE,
                                         ...)

  ndim <- length(lowerLimit)

  if (is.null(which)) {
    fnIndices <- seq_along(fns)
  } else {
    fnIndices <- na.omit(match(which, names(fns)))
  }
  fnList <- lapply(names(fns)[fnIndices], function(x) call(x))

  argList <- c(fnList, iterations = iterations, check = FALSE)
  result <- do.call(bench::mark, args = argList)
  d <- result[seq_along(fnIndices), 1:9]
  d$expression <- fns[fnIndices]
  d
}
```

We reel off the timing runs.

## Example 1

``` r
func <- function(x) sin(x[1]) * cos(x[2]) * exp(x[3])
func.v <- function(x) {
    matrix(apply(x, 2, function(z) sin(z[1]) * cos(z[2]) * exp(z[3])), ncol = ncol(x))
}

d <- harness(f = func, fv = func.v,
             lowerLimit = rep(0, 3),
             upperLimit = rep(1, 3),
             tol = 1e-5,
             iterations = 100)[, 1:9]

knitr::kable(d, digits = 3, row.names = FALSE)
```

| expression                     |      min |   median |  itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|---------:|---------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 274.04µs | 284.03µs | 3449.682 |    64.4KB | 34.845 |    99 |    1 |     28.7ms |
| Vectorized Hcubature           | 392.46µs |  420.8µs | 2325.907 |    54.5KB |  0.000 |   100 |    0 |       43ms |
| Non-vectorized Pcubature       | 880.73µs | 899.65µs | 1092.113 |    39.2KB | 33.777 |    97 |    3 |     88.8ms |
| Vectorized Pcubature           |   1.18ms |   1.21ms |  817.168 |    58.4KB | 16.677 |    98 |    2 |    119.9ms |
| Non-vectorized cubature::cuhre |  581.9µs | 607.26µs | 1633.453 |    40.4KB | 33.336 |    98 |    2 |       60ms |
| Vectorized cubature::cuhre     | 596.74µs | 627.18µs | 1587.197 |    39.8KB | 16.032 |    99 |    1 |     62.4ms |

## Multivariate normal

Using `cubature`, we evaluate $$\int_{R}\phi(x)dx$$ where $\phi(x)$ is
the three-dimensional multivariate normal density with mean 0, and
variance $$\Sigma = \left( \begin{array}{rrr}
1 & \frac{3}{5} & \frac{1}{3} \\
\frac{3}{5} & 1 & \frac{11}{15} \\
\frac{1}{3} & \frac{11}{15} & 1
\end{array} \right)$$ and $R$ is
$\left\lbrack - \frac{1}{2},1 \right\rbrack \times \left\lbrack - \frac{1}{2},4 \right\rbrack \times \left\lbrack - \frac{1}{2},2 \right\rbrack.$

We construct a scalar function (`my_dmvnorm`) and a vector analogue
(`my_dmvnorm_v`). First the functions.

``` r
m <- 3
sigma <- diag(3)
sigma[2,1] <- sigma[1, 2] <- 3/5 ; sigma[3,1] <- sigma[1, 3] <- 1/3
sigma[3,2] <- sigma[2, 3] <- 11/15
logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
my_dmvnorm <- function (x, mean, sigma, logdet) {
    x <- matrix(x, ncol = length(x))
    distval <- stats::mahalanobis(x, center = mean, cov = sigma)
    exp(-(3 * log(2 * pi) + logdet + distval)/2)
}

my_dmvnorm_v <- function (x, mean, sigma, logdet) {
    distval <- stats::mahalanobis(t(x), center = mean, cov = sigma)
    exp(matrix(-(3 * log(2 * pi) + logdet + distval)/2, ncol = ncol(x)))
}
```

Now the timing.

``` r
d <- harness(f = my_dmvnorm, fv = my_dmvnorm_v,
             lowerLimit = rep(-0.5, 3),
             upperLimit = c(1, 4, 2),
             tol = 1e-5,
             iterations = 10,
             mean = rep(0, m), sigma = sigma, logdet = logdet)
knitr::kable(d, digits = 3)
```

| expression                     |      min |   median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|---------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 830.49ms | 862.21ms |   1.157 |  139.68KB | 18.403 |    10 |  159 |      8.64s |
| Vectorized Hcubature           |   1.78ms |   2.09ms | 442.878 |    1.89MB |  0.000 |    10 |    0 |    22.58ms |
| Non-vectorized Pcubature       | 364.48ms | 366.37ms |   2.674 |        0B | 18.186 |    10 |   68 |      3.74s |
| Vectorized Pcubature           |   1.28ms |   1.32ms | 734.940 |  810.26KB |  0.000 |    10 |    0 |    13.61ms |
| Non-vectorized cubature::cuhre |  349.5ms | 351.83ms |   2.841 |        0B | 18.470 |    10 |   65 |      3.52s |
| Vectorized cubature::cuhre     |    3.3ms |    3.4ms | 294.652 |  898.41KB |  0.000 |    10 |    0 |    33.94ms |

The effect of vectorization is huge. So it makes sense for users to
vectorize their integrands as much as possible for efficiency.

Furthermore, for this particular example, we expect
[`mvtnorm::pmvnorm`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html) to do
pretty well since it is specialized for the multivariate normal. The
vectorized versions of `hcubature` and `pcubature` seem competitive and
in some cases better, as the table below shows.

``` r
library(mvtnorm)
g1 <- function() pmvnorm(lower = rep(-0.5, m),
                                  upper = c(1, 4, 2), mean = rep(0, m), corr = sigma,
                                  alg = Miwa(), abseps = 1e-5, releps = 1e-5)
g2 <- function() pmvnorm(lower = rep(-0.5, m),
                         upper = c(1, 4, 2), mean = rep(0, m), corr = sigma,
                         alg = GenzBretz(), abseps = 1e-5, releps = 1e-5)
g3 <- function() pmvnorm(lower = rep(-0.5, m),
                         upper = c(1, 4, 2), mean = rep(0, m), corr = sigma,
                         alg = TVPACK(), abseps = 1e-5, releps = 1e-5)

knitr::kable(bench::mark(g1(), g2(), g3(), iterations = 20, check = FALSE)[, 1:9],
             digits = 3, row.names = FALSE)
```

| expression |   min | median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-----------|------:|-------:|--------:|----------:|-------:|------:|-----:|-----------:|
| g1()       | 775µs | 1.38ms | 731.408 |  355.13KB |      0 |    20 |    0 |     27.3ms |
| g2()       | 751µs | 1.38ms | 750.304 |    2.49KB |      0 |    20 |    0 |     26.7ms |
| g3()       | 755µs | 1.39ms | 725.739 |    2.49KB |      0 |    20 |    0 |     27.6ms |

## Product of cosines

``` r
testFn0 <- function(x) prod(cos(x))
testFn0_v <- function(x) matrix(apply(x, 2, function(z) prod(cos(z))), ncol = ncol(x))

d <- harness(f = testFn0, fv = testFn0_v,
             lowerLimit = rep(0, 2), upperLimit = rep(1, 2), iterations = 1000)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  37.6µs |  39.5µs | 24300.767 |    7.66KB | 48.699 |   998 |    2 |     41.1ms |
| Vectorized Hcubature           |  66.7µs |  69.7µs | 13826.205 |    1.16KB | 13.840 |   999 |    1 |     72.3ms |
| Non-vectorized Pcubature       |  51.5µs |  53.7µs | 17873.315 |        0B | 35.818 |   998 |    2 |     55.8ms |
| Vectorized Pcubature           | 118.3µs | 124.6µs |  7776.282 |   18.68KB | 23.399 |   997 |    3 |    128.2ms |
| Non-vectorized cubature::cuhre | 329.6µs | 345.1µs |  2877.861 |        0B | 26.136 |   991 |    9 |    344.4ms |
| Vectorized cubature::cuhre     | 365.4µs | 383.7µs |  2571.760 |   16.38KB | 25.977 |   990 |   10 |    384.9ms |

## Gaussian function

``` r
testFn1 <- function(x) {
    val <- sum(((1 - x) / x)^2)
    scale <- prod((2 / sqrt(pi)) / x^2)
    exp(-val) * scale
}

testFn1_v <- function(x) {
    val <- matrix(apply(x, 2, function(z) sum(((1 - z) / z)^2)), ncol(x))
    scale <- matrix(apply(x, 2, function(z) prod((2 / sqrt(pi)) / z^2)), ncol(x))
    exp(-val) * scale
}

d <- harness(f = testFn1, fv = testFn1_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 10)

knitr::kable(d, digits = 3)
```

| expression                     |      min |  median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|--------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |   2.96ms |  3.05ms |   326.231 |    67.3KB | 36.248 |     9 |    1 |     27.6ms |
| Vectorized Hcubature           |   4.92ms |  4.99ms |   199.126 |   290.5KB | 22.125 |     9 |    1 |     45.2ms |
| Non-vectorized Pcubature       |  78.37µs | 80.27µs | 11573.811 |        0B |  0.000 |    10 |    0 |      864µs |
| Vectorized Pcubature           | 157.32µs |   162µs |  5892.604 |     4.1KB |  0.000 |    10 |    0 |      1.7ms |
| Non-vectorized cubature::cuhre |  13.72ms | 13.91ms |    71.650 |        0B | 47.767 |     6 |    4 |     83.7ms |
| Vectorized cubature::cuhre     |   20.4ms | 20.45ms |    48.800 |   971.5KB | 48.800 |     5 |    5 |    102.5ms |

## Discontinuous function

``` r
testFn2 <- function(x) {
    radius <- 0.50124145262344534123412
    ifelse(sum(x * x) < radius * radius, 1, 0)
}

testFn2_v <- function(x) {
    radius <- 0.50124145262344534123412
    matrix(apply(x, 2, function(z) ifelse(sum(z * z) < radius * radius, 1, 0)), ncol = ncol(x))
}

d <- harness(which = c("hc", "hc.v", "cc", "cc_v"),
             f = testFn2, fv = testFn2_v,
             lowerLimit = rep(0, 2), upperLimit = rep(1, 2), iterations = 10)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  42.4ms |  44.5ms |  18.811 |    17.7KB | 30.098 |    10 |   16 |    531.6ms |
| Vectorized Hcubature           |  49.2ms |  49.7ms |  19.855 |  1011.8KB | 25.812 |    10 |   13 |   503.64ms |
| Non-vectorized cubature::cuhre | 823.2ms | 832.6ms |   1.161 |        0B | 28.221 |    10 |  243 |      8.61s |
| Vectorized cubature::cuhre     | 905.6ms | 920.8ms |   1.076 |    21.2MB | 25.275 |    10 |  235 |       9.3s |

## A simple polynomial (product of coordinates)

``` r
testFn3 <- function(x) prod(2 * x)
testFn3_v <- function(x) matrix(apply(x, 2, function(z) prod(2 * z)), ncol = ncol(x))

d <- harness(f = testFn3, fv = testFn3_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 50)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  58.8µs |    62µs | 15377.525 |    6.75KB |  0.000 |    50 |    0 |     3.25ms |
| Vectorized Hcubature           |  90.5µs |  95.5µs |  9999.695 |    2.91KB |  0.000 |    50 |    0 |        5ms |
| Non-vectorized Pcubature       |  53.5µs |  55.1µs | 17306.801 |        0B |  0.000 |    50 |    0 |     2.89ms |
| Vectorized Pcubature           |  85.5µs |  88.3µs | 10809.590 |   19.12KB |  0.000 |    50 |    0 |     4.63ms |
| Non-vectorized cubature::cuhre | 634.1µs | 661.5µs |  1498.883 |        0B | 30.589 |    49 |    1 |    32.69ms |
| Vectorized cubature::cuhre     | 637.2µs | 663.4µs |  1488.937 |   39.84KB | 30.386 |    49 |    1 |    32.91ms |

## Gaussian centered at $\frac{1}{2}$

``` r
testFn4 <- function(x) {
    a <- 0.1
    s <- sum((x - 0.5)^2)
    ((2 / sqrt(pi)) / (2. * a))^length(x) * exp (-s / (a * a))
}

testFn4_v <- function(x) {
    a <- 0.1
    r <- apply(x, 2, function(z) {
        s <- sum((z - 0.5)^2)
        ((2 / sqrt(pi)) / (2. * a))^length(z) * exp (-s / (a * a))
    })
    matrix(r, ncol = ncol(x))
}

d <- harness(f = testFn4, fv = testFn4_v,
             lowerLimit = rep(0, 2), upperLimit = rep(1, 2), iterations = 20)
knitr::kable(d, digits = 3)
```

| expression                     |    min | median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|-------:|-------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 1.35ms | 1.39ms | 716.454 |    85.2KB | 37.708 |    19 |    1 |     26.5ms |
| Vectorized Hcubature           | 1.76ms | 1.82ms | 544.820 |   147.2KB | 28.675 |    19 |    1 |     34.9ms |
| Non-vectorized Pcubature       | 2.04ms | 2.09ms | 477.496 |        0B | 25.131 |    19 |    1 |     39.8ms |
| Vectorized Pcubature           | 2.57ms | 2.63ms | 379.070 |    68.9KB | 19.951 |    19 |    1 |     50.1ms |
| Non-vectorized cubature::cuhre | 3.21ms | 3.26ms | 303.982 |        0B | 33.776 |    18 |    2 |     59.2ms |
| Vectorized cubature::cuhre     | 3.72ms | 3.79ms | 247.447 |   125.6KB | 27.494 |    18 |    2 |     72.7ms |

## Double Gaussian

``` r
testFn5 <- function(x) {
    a <- 0.1
    s1 <- sum((x - 1 / 3)^2)
    s2 <- sum((x - 2 / 3)^2)
    0.5 * ((2 / sqrt(pi)) / (2. * a))^length(x) * (exp(-s1 / (a * a)) + exp(-s2 / (a * a)))
}
testFn5_v <- function(x) {
    a <- 0.1
    r <- apply(x, 2, function(z) {
        s1 <- sum((z - 1 / 3)^2)
        s2 <- sum((z - 2 / 3)^2)
        0.5 * ((2 / sqrt(pi)) / (2. * a))^length(z) * (exp(-s1 / (a * a)) + exp(-s2 / (a * a)))
    })
    matrix(r, ncol = ncol(x))
}

d <- harness(f = testFn5, fv = testFn5_v,
             lowerLimit = rep(0, 2), upperLimit = rep(1, 2), iterations = 20)
knitr::kable(d, digits = 3)
```

| expression                     |    min | median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|-------:|-------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 3.66ms | 3.75ms | 262.533 |     133KB | 29.170 |    18 |    2 |     68.6ms |
| Vectorized Hcubature           |  4.7ms |  4.8ms | 194.556 |     249KB | 34.333 |    17 |    3 |     87.4ms |
| Non-vectorized Pcubature       | 2.54ms |  2.6ms | 319.219 |        0B | 16.801 |    19 |    1 |     59.5ms |
| Vectorized Pcubature           |  3.3ms | 3.42ms | 292.346 |      69KB | 32.483 |    18 |    2 |     61.6ms |
| Non-vectorized cubature::cuhre |    7ms | 7.15ms | 139.260 |        0B | 34.815 |    16 |    4 |    114.9ms |
| Vectorized cubature::cuhre     | 8.38ms | 8.51ms | 117.014 |     224KB | 29.253 |    16 |    4 |    136.7ms |

## Tsuda’s example

``` r
testFn6 <- function(x) {
    a <- (1 + sqrt(10.0)) / 9.0
    prod( a / (a + 1) * ((a + 1) / (a + x))^2)
}

testFn6_v <- function(x) {
    a <- (1 + sqrt(10.0)) / 9.0
    r <- apply(x, 2, function(z) prod( a / (a + 1) * ((a + 1) / (a + z))^2))
    matrix(r, ncol = ncol(x))
}

d <- harness(f = testFn6, fv = testFn6_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 20)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  1.57ms |  1.61ms | 617.600 |    64.6KB | 32.505 |    19 |    1 |     30.8ms |
| Vectorized Hcubature           |  2.14ms |  2.19ms | 457.482 |     156KB | 24.078 |    19 |    1 |     41.5ms |
| Non-vectorized Pcubature       |  8.27ms |  8.37ms | 118.816 |        0B | 39.605 |    15 |    5 |    126.2ms |
| Vectorized Pcubature           | 10.25ms | 10.41ms |  96.128 |   386.2KB | 32.043 |    15 |    5 |      156ms |
| Non-vectorized cubature::cuhre |  4.33ms |  4.42ms | 226.815 |        0B | 25.202 |    18 |    2 |     79.4ms |
| Vectorized cubature::cuhre     |  4.72ms |  4.83ms | 205.075 |   225.8KB | 36.190 |    17 |    3 |     82.9ms |

## Morokoff & Caflisch example

``` r
testFn7 <- function(x) {
    n <- length(x)
    p <- 1/n
    (1 + p)^n * prod(x^p)
}
testFn7_v <- function(x) {
    matrix(apply(x, 2, function(z) {
        n <- length(z)
        p <- 1/n
        (1 + p)^n * prod(z^p)
    }), ncol = ncol(x))
}

d <- harness(f = testFn7, fv = testFn7_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 20)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |   3.4ms |  3.46ms | 264.467 |   32.89KB | 26.447 |    20 |    2 |     75.6ms |
| Vectorized Hcubature           |  4.04ms |   4.2ms | 223.250 |  205.61KB | 22.325 |    20 |    2 |     89.6ms |
| Non-vectorized Pcubature       |  8.36ms |  8.59ms | 103.997 |        0B | 31.199 |    20 |    6 |    192.3ms |
| Vectorized Pcubature           |  9.59ms |  9.83ms |  94.387 |  386.24KB | 23.597 |    20 |    5 |    211.9ms |
| Non-vectorized cubature::cuhre | 44.05ms | 44.64ms |  22.207 |        0B | 25.538 |    20 |   23 |    900.6ms |
| Vectorized cubature::cuhre     | 44.01ms | 45.53ms |  21.890 |    2.04MB | 24.079 |    20 |   22 |    913.6ms |

## Wang-Landau sampling 1d, 2d examples

``` r
I.1d <- function(x) {
    sin(4 * x) *
        x * ((x * ( x * (x * x - 4) + 1) - 1))
}
I.1d_v <- function(x) {
    matrix(apply(x, 2, function(z)
        sin(4 * z) *
        z * ((z * ( z * (z * z - 4) + 1) - 1))),
        ncol = ncol(x))
}
d <- harness(f = I.1d, fv = I.1d_v,
             lowerLimit = -2, upperLimit = 2, iterations = 100)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 127.9µs |   135µs |  6508.061 |    53.7KB | 65.738 |    99 |    1 |    15.21ms |
| Vectorized Hcubature           | 208.7µs | 216.2µs |  4486.378 |    68.5KB |  0.000 |   100 |    0 |    22.29ms |
| Non-vectorized Pcubature       |  55.2µs |  57.2µs | 16762.951 |        0B |  0.000 |   100 |    0 |     5.97ms |
| Vectorized Pcubature           | 156.2µs | 164.5µs |  5776.581 |        0B | 58.349 |    99 |    1 |    17.14ms |
| Non-vectorized cubature::cuhre | 220.3µs | 226.6µs |  4348.902 |        0B |  0.000 |   100 |    0 |    22.99ms |
| Vectorized cubature::cuhre     | 495.9µs | 524.2µs |  1886.810 |        0B | 38.506 |    98 |    2 |    51.94ms |

``` r
I.2d <- function(x) {
    x1 <- x[1]; x2 <- x[2]
    sin(4 * x1 + 1) * cos(4 * x2) * x1 * (x1 * (x1 * x1)^2 - x2 * (x2 * x2 - x1) +2)
}
I.2d_v <- function(x) {
    matrix(apply(x, 2,
                 function(z) {
                     x1 <- z[1]; x2 <- z[2]
                     sin(4 * x1 + 1) * cos(4 * x2) * x1 * (x1 * (x1 * x1)^2 - x2 * (x2 * x2 - x1) +2)
                 }),
           ncol = ncol(x))
}
d <- harness(f = I.2d, fv = I.2d_v,
             lowerLimit = rep(-1, 2), upperLimit = rep(1, 2), iterations = 100)
knitr::kable(d, digits = 3)
```

| expression                     |      min |   median |  itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|---------:|---------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |   4.47ms |   4.55ms |  218.655 |    78.4KB | 35.595 |    86 |   14 |    393.3ms |
| Vectorized Hcubature           |   5.33ms |   5.43ms |  183.737 |   304.7KB | 29.911 |    86 |   14 |    468.1ms |
| Non-vectorized Pcubature       | 414.57µs | 435.87µs | 2287.366 |        0B | 23.105 |    99 |    1 |     43.3ms |
| Vectorized Pcubature           | 602.59µs | 630.69µs | 1573.608 |    18.3KB | 15.895 |    99 |    1 |     62.9ms |
| Non-vectorized cubature::cuhre |   1.22ms |   1.25ms |  796.377 |        0B | 24.630 |    97 |    3 |    121.8ms |
| Vectorized cubature::cuhre     |   1.32ms |   1.36ms |  730.142 |    60.1KB | 30.423 |    96 |    4 |    131.5ms |

## Summary

The recommendations are:

1.  **Vectorize your integrand.** The time spent in so doing pays back
    enormously on all but the most trivial problems. This is easy to do
    and the examples above show how.

2.  **Vectorized `hcubature` is a good default starting point.** It
    handles the widest range of problems well.

3.  **For smooth integrands in low dimensions ($\leq 3$), try
    `pcubature`.** It can be substantially faster than `hcubature` when
    its assumptions hold. Experiment before committing in a production
    setting, and be aware of the sampling-density cliff documented in
    the Get Started vignette — when in doubt, cross-check with `cuhre`.

4.  **For non-smooth or localized integrands, prefer
    `cubintegrate(method = "cuhre")`** or at least
    `hcubature(..., robust = TRUE)`. See the Robustness section of the
    [Get Started](cubature.md) vignette.

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
    ## [1] mvtnorm_1.3-6  cubature_2.2.0 bench_1.1.4   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.7.3       cli_3.6.6         knitr_1.51        rlang_1.2.0      
    ##  [5] xfun_0.57         textshaping_1.0.5 jsonlite_2.0.0    glue_1.8.0       
    ##  [9] htmltools_0.5.9   ragg_1.5.2        sass_0.4.10       rmarkdown_2.31   
    ## [13] tibble_3.3.1      evaluate_1.0.5    jquerylib_0.1.4   fastmap_1.2.0    
    ## [17] profmem_0.7.0     yaml_2.3.12       lifecycle_1.0.5   compiler_4.5.3   
    ## [21] fs_2.0.1          pkgconfig_2.0.3   Rcpp_1.1.1        systemfonts_1.3.2
    ## [25] digest_0.6.39     R6_2.6.1          pillar_1.11.1     magrittr_2.0.5   
    ## [29] bslib_0.10.0      tools_4.5.3       pkgdown_2.2.0     cachem_1.1.0     
    ## [33] desc_1.4.3
