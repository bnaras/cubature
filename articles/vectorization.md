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
| Non-vectorized Hcubature       | 273.23µs | 286.79µs | 3460.028 |    64.4KB | 34.950 |    99 |    1 |     28.6ms |
| Vectorized Hcubature           | 384.76µs | 411.39µs | 2379.679 |    54.5KB |  0.000 |   100 |    0 |       42ms |
| Non-vectorized Pcubature       | 865.24µs |  885.7µs | 1112.672 |    39.2KB | 34.413 |    97 |    3 |     87.2ms |
| Vectorized Pcubature           |   1.16ms |   1.19ms |  833.805 |    58.4KB | 17.016 |    98 |    2 |    117.5ms |
| Non-vectorized cubature::cuhre | 579.87µs | 601.28µs | 1646.726 |    40.4KB | 33.607 |    98 |    2 |     59.5ms |
| Vectorized cubature::cuhre     | 582.62µs | 615.66µs | 1612.229 |    39.8KB | 16.285 |    99 |    1 |     61.4ms |

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
| Non-vectorized Hcubature       | 816.85ms | 820.59ms |   1.200 |  139.68KB | 19.084 |    10 |  159 |      8.33s |
| Vectorized Hcubature           |   1.79ms |   2.15ms | 448.431 |    1.89MB |  0.000 |    10 |    0 |     22.3ms |
| Non-vectorized Pcubature       |  352.7ms | 354.37ms |   2.760 |        0B | 18.769 |    10 |   68 |      3.62s |
| Vectorized Pcubature           |   1.26ms |   1.28ms | 762.980 |  810.26KB |  0.000 |    10 |    0 |    13.11ms |
| Non-vectorized cubature::cuhre | 338.58ms | 341.54ms |   2.927 |        0B | 19.026 |    10 |   65 |      3.42s |
| Vectorized cubature::cuhre     |   3.25ms |   3.28ms | 302.530 |  898.41KB |  0.000 |    10 |    0 |    33.05ms |

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
| g1()       | 757µs | 1.38ms | 734.166 |  355.13KB |      0 |    20 |    0 |     27.2ms |
| g2()       | 765µs | 1.38ms | 758.771 |    2.49KB |      0 |    20 |    0 |     26.4ms |
| g3()       | 752µs | 1.38ms | 739.555 |    2.49KB |      0 |    20 |    0 |       27ms |

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
| Non-vectorized Hcubature       |  37.1µs |    39µs | 24574.794 |    7.66KB | 49.248 |   998 |    2 |     40.6ms |
| Vectorized Hcubature           |  66.3µs |  69.2µs | 13976.476 |    1.16KB | 13.990 |   999 |    1 |     71.5ms |
| Non-vectorized Pcubature       |  50.6µs |  52.7µs | 18323.884 |        0B | 36.721 |   998 |    2 |     54.5ms |
| Vectorized Pcubature           | 116.6µs | 122.1µs |  7884.054 |   18.68KB | 23.723 |   997 |    3 |    126.5ms |
| Non-vectorized cubature::cuhre | 328.8µs | 339.6µs |  2905.470 |        0B | 26.387 |   991 |    9 |    341.1ms |
| Vectorized cubature::cuhre     | 360.8µs | 375.2µs |  2634.434 |   16.38KB | 26.610 |   990 |   10 |    375.8ms |

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

| expression                     |      min |   median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|---------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |   2.93ms |   2.97ms |   333.770 |    67.3KB | 37.086 |     9 |    1 |    26.96ms |
| Vectorized Hcubature           |   4.86ms |      5ms |   200.125 |   290.5KB | 22.236 |     9 |    1 |    44.97ms |
| Non-vectorized Pcubature       |  78.06µs |  82.68µs | 11230.483 |        0B |  0.000 |    10 |    0 |   890.43µs |
| Vectorized Pcubature           | 155.52µs | 161.27µs |  6031.830 |     4.1KB |  0.000 |    10 |    0 |     1.66ms |
| Non-vectorized cubature::cuhre |   13.6ms |  13.76ms |    72.472 |        0B | 48.315 |     6 |    4 |    82.79ms |
| Vectorized cubature::cuhre     |  20.17ms |  20.28ms |    49.288 |   971.5KB | 49.288 |     5 |    5 |   101.44ms |

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
| Non-vectorized Hcubature       |  40.9ms |  42.6ms |  19.127 |    17.7KB | 30.603 |    10 |   16 |   522.83ms |
| Vectorized Hcubature           |  47.8ms |  48.4ms |  20.566 |  1011.8KB | 26.736 |    10 |   13 |   486.24ms |
| Non-vectorized cubature::cuhre | 789.4ms | 798.3ms |   1.235 |        0B | 30.130 |    10 |  244 |       8.1s |
| Vectorized cubature::cuhre     | 890.5ms | 901.4ms |   1.102 |    21.2MB | 25.794 |    10 |  234 |      9.07s |

## A simple polynomial (product of coordinates)

``` r
testFn3 <- function(x) prod(2 * x)
testFn3_v <- function(x) matrix(apply(x, 2, function(z) prod(2 * z)), ncol = ncol(x))

d <- harness(f = testFn3, fv = testFn3_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 50)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median |   itr/sec | mem_alloc |  gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|----------:|----------:|--------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  59.5µs |  61.4µs | 15355.364 |    6.75KB |   0.000 |    50 |    0 |     3.26ms |
| Vectorized Hcubature           |  92.5µs |  96.8µs | 10040.206 |    2.91KB |   0.000 |    50 |    0 |     4.98ms |
| Non-vectorized Pcubature       |  55.6µs |  58.3µs | 15991.049 |        0B | 326.348 |    49 |    1 |     3.06ms |
| Vectorized Pcubature           |  84.4µs |  88.2µs | 10918.069 |   19.12KB |   0.000 |    50 |    0 |     4.58ms |
| Non-vectorized cubature::cuhre | 637.2µs |   663µs |  1500.595 |        0B |  30.624 |    49 |    1 |    32.65ms |
| Vectorized cubature::cuhre     | 637.8µs | 664.2µs |  1500.822 |   39.84KB |   0.000 |    50 |    0 |    33.31ms |

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
| Non-vectorized Hcubature       | 1.35ms | 1.39ms | 718.005 |    85.2KB | 37.790 |    19 |    1 |     26.5ms |
| Vectorized Hcubature           | 1.78ms | 1.85ms | 541.106 |   147.2KB | 28.479 |    19 |    1 |     35.1ms |
| Non-vectorized Pcubature       | 2.05ms | 2.08ms | 479.438 |        0B | 25.234 |    19 |    1 |     39.6ms |
| Vectorized Pcubature           | 2.56ms | 2.62ms | 380.795 |    68.9KB | 20.042 |    19 |    1 |     49.9ms |
| Non-vectorized cubature::cuhre | 3.21ms | 3.24ms | 307.896 |        0B | 34.211 |    18 |    2 |     58.5ms |
| Vectorized cubature::cuhre     | 3.75ms | 3.81ms | 262.196 |   125.6KB | 29.133 |    18 |    2 |     68.7ms |

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
| Non-vectorized Hcubature       | 3.61ms | 3.68ms | 271.656 |     133KB | 30.184 |    18 |    2 |     66.3ms |
| Vectorized Hcubature           | 4.56ms | 4.64ms | 214.092 |     249KB | 37.781 |    17 |    3 |     79.4ms |
| Non-vectorized Pcubature       |  2.5ms | 2.53ms | 393.154 |        0B | 20.692 |    19 |    1 |     48.3ms |
| Vectorized Pcubature           | 3.24ms | 3.31ms | 298.333 |      69KB | 33.148 |    18 |    2 |     60.3ms |
| Non-vectorized cubature::cuhre | 6.97ms | 7.05ms | 141.376 |        0B | 35.344 |    16 |    4 |    113.2ms |
| Vectorized cubature::cuhre     | 8.32ms | 8.43ms | 116.048 |     224KB | 29.012 |    16 |    4 |    137.9ms |

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
| Non-vectorized Hcubature       |  1.56ms |   1.6ms | 624.762 |    64.6KB | 32.882 |    19 |    1 |     30.4ms |
| Vectorized Hcubature           |  2.11ms |  2.14ms | 464.839 |     156KB | 24.465 |    19 |    1 |     40.9ms |
| Non-vectorized Pcubature       |  8.12ms |   8.2ms | 121.136 |        0B | 40.379 |    15 |    5 |    123.8ms |
| Vectorized Pcubature           | 10.21ms | 10.38ms |  96.451 |   386.2KB | 32.150 |    15 |    5 |    155.5ms |
| Non-vectorized cubature::cuhre |  4.31ms |  4.38ms | 227.697 |        0B | 25.300 |    18 |    2 |     79.1ms |
| Vectorized cubature::cuhre     |  4.73ms |  4.82ms | 207.175 |   225.8KB | 36.560 |    17 |    3 |     82.1ms |

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
| Non-vectorized Hcubature       |  3.34ms |  3.41ms | 275.491 |   32.89KB | 27.549 |    20 |    2 |     72.6ms |
| Vectorized Hcubature           |  3.99ms |  4.08ms | 233.348 |  205.61KB | 23.335 |    20 |    2 |     85.7ms |
| Non-vectorized Pcubature       |  8.25ms |  8.45ms | 108.773 |        0B | 32.632 |    20 |    6 |    183.9ms |
| Vectorized Pcubature           |  9.39ms |  9.63ms |  98.014 |  386.24KB | 24.503 |    20 |    5 |    204.1ms |
| Non-vectorized cubature::cuhre | 42.55ms | 42.85ms |  23.085 |        0B | 26.548 |    20 |   23 |    866.4ms |
| Vectorized cubature::cuhre     | 42.62ms | 43.22ms |  22.970 |    2.04MB | 25.267 |    20 |   22 |    870.7ms |

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

| expression                     |   min |  median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|------:|--------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 130µs | 133.3µs |  7320.174 |    53.7KB |  0.000 |   100 |    0 |    13.66ms |
| Vectorized Hcubature           | 213µs | 220.9µs |  4371.734 |    68.5KB | 44.159 |    99 |    1 |    22.64ms |
| Non-vectorized Pcubature       |  55µs |  57.1µs | 16713.644 |        0B |  0.000 |   100 |    0 |     5.98ms |
| Vectorized Pcubature           | 157µs | 163.1µs |  5947.573 |        0B |  0.000 |   100 |    0 |    16.81ms |
| Non-vectorized cubature::cuhre | 226µs | 232.5µs |  4184.584 |        0B | 42.269 |    99 |    1 |    23.66ms |
| Vectorized cubature::cuhre     | 495µs | 516.4µs |  1927.086 |        0B | 19.466 |    99 |    1 |    51.37ms |

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
| Non-vectorized Hcubature       |   4.51ms |    4.6ms |  217.273 |    78.4KB | 35.370 |    86 |   14 |    395.8ms |
| Vectorized Hcubature           |   5.32ms |    5.4ms |  184.589 |   304.7KB | 27.582 |    87 |   13 |    471.3ms |
| Non-vectorized Pcubature       | 413.79µs | 429.15µs | 2318.957 |        0B | 23.424 |    99 |    1 |     42.7ms |
| Vectorized Pcubature           | 595.41µs | 624.34µs | 1582.770 |    18.3KB | 32.301 |    98 |    2 |     61.9ms |
| Non-vectorized cubature::cuhre |   1.22ms |   1.24ms |  804.141 |        0B | 24.870 |    97 |    3 |    120.6ms |
| Vectorized cubature::cuhre     |   1.32ms |   1.35ms |  733.460 |    60.1KB | 22.684 |    97 |    3 |    132.2ms |

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
