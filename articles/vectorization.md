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
| Non-vectorized Hcubature       | 272.22µs |  282.3µs | 3479.881 |    64.4KB | 35.150 |    99 |    1 |     28.4ms |
| Vectorized Hcubature           | 387.94µs | 416.42µs | 2371.814 |    54.5KB |  0.000 |   100 |    0 |     42.2ms |
| Non-vectorized Pcubature       | 853.05µs | 881.33µs | 1119.147 |    39.2KB | 34.613 |    97 |    3 |     86.7ms |
| Vectorized Pcubature           |   1.16ms |   1.19ms |  833.351 |    58.4KB | 17.007 |    98 |    2 |    117.6ms |
| Non-vectorized cubature::cuhre | 576.69µs | 597.38µs | 1660.635 |    40.4KB | 33.891 |    98 |    2 |       59ms |
| Vectorized cubature::cuhre     | 577.11µs | 612.79µs | 1626.880 |    39.8KB | 16.433 |    99 |    1 |     60.9ms |

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
| Non-vectorized Hcubature       | 815.37ms | 824.87ms |   1.195 |  139.68KB | 19.002 |    10 |  159 |      8.37s |
| Vectorized Hcubature           |   1.79ms |   2.12ms | 443.194 |    1.89MB |  0.000 |    10 |    0 |    22.56ms |
| Non-vectorized Pcubature       | 351.68ms |  354.9ms |   2.771 |        0B | 18.845 |    10 |   68 |      3.61s |
| Vectorized Pcubature           |   1.26ms |    1.3ms | 761.638 |  810.26KB |  0.000 |    10 |    0 |    13.13ms |
| Non-vectorized cubature::cuhre | 338.48ms | 341.42ms |   2.932 |        0B | 19.056 |    10 |   65 |      3.41s |
| Vectorized cubature::cuhre     |   3.23ms |   3.28ms | 301.551 |  898.41KB |  0.000 |    10 |    0 |    33.16ms |

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
| g1()       | 764µs | 1.38ms | 733.532 |  355.13KB |      0 |    20 |    0 |     27.3ms |
| g2()       | 766µs | 1.37ms | 755.885 |    2.49KB |      0 |    20 |    0 |     26.5ms |
| g3()       | 751µs | 1.38ms | 737.770 |    2.49KB |      0 |    20 |    0 |     27.1ms |

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
| Non-vectorized Hcubature       |  37.7µs |  39.4µs | 24468.103 |    7.66KB | 49.034 |   998 |    2 |     40.8ms |
| Vectorized Hcubature           |  65.6µs |  69.1µs | 13995.017 |    1.16KB | 14.009 |   999 |    1 |     71.4ms |
| Non-vectorized Pcubature       |  50.6µs |  52.8µs | 18215.148 |        0B | 36.503 |   998 |    2 |     54.8ms |
| Vectorized Pcubature           | 118.1µs | 124.1µs |  7828.974 |   18.68KB | 23.558 |   997 |    3 |    127.3ms |
| Non-vectorized cubature::cuhre | 329.5µs | 341.6µs |  2901.090 |        0B | 26.347 |   991 |    9 |    341.6ms |
| Vectorized cubature::cuhre     |   358µs | 375.5µs |  2633.551 |   16.38KB | 26.602 |   990 |   10 |    375.9ms |

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
| Non-vectorized Hcubature       |   2.88ms |   2.93ms |   339.354 |    67.3KB | 37.706 |     9 |    1 |    26.52ms |
| Vectorized Hcubature           |   4.81ms |   4.93ms |   203.249 |   290.5KB | 22.583 |     9 |    1 |    44.28ms |
| Non-vectorized Pcubature       |  77.99µs |  79.56µs | 11659.037 |        0B |  0.000 |    10 |    0 |    857.7µs |
| Vectorized Pcubature           | 154.47µs | 161.69µs |  5944.676 |     4.1KB |  0.000 |    10 |    0 |     1.68ms |
| Non-vectorized cubature::cuhre |  13.49ms |  13.78ms |    72.725 |        0B | 48.483 |     6 |    4 |     82.5ms |
| Vectorized cubature::cuhre     |  20.26ms |  20.42ms |    49.006 |   971.5KB | 49.006 |     5 |    5 |   102.03ms |

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
| Non-vectorized Hcubature       |  40.9ms |  42.3ms |  19.263 |    17.7KB | 30.821 |    10 |   16 |   519.12ms |
| Vectorized Hcubature           |  46.8ms |  48.6ms |  20.681 |  1011.8KB | 26.885 |    10 |   13 |   483.54ms |
| Non-vectorized cubature::cuhre | 795.4ms | 808.3ms |   1.216 |        0B | 29.681 |    10 |  244 |      8.22s |
| Vectorized cubature::cuhre     | 884.9ms | 891.9ms |   1.112 |    21.2MB | 26.022 |    10 |  234 |      8.99s |

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
| Non-vectorized Hcubature       |  59.1µs |  60.9µs | 15906.021 |    6.75KB |   0.000 |    50 |    0 |     3.14ms |
| Vectorized Hcubature           |  92.4µs |   100µs |  9494.738 |    2.91KB | 193.770 |    49 |    1 |     5.16ms |
| Non-vectorized Pcubature       |  53.5µs |  54.8µs | 17490.903 |        0B |   0.000 |    50 |    0 |     2.86ms |
| Vectorized Pcubature           |  84.1µs |  87.8µs | 11012.801 |   19.12KB |   0.000 |    50 |    0 |     4.54ms |
| Non-vectorized cubature::cuhre | 628.4µs | 662.9µs |  1499.325 |        0B |  30.598 |    49 |    1 |    32.68ms |
| Vectorized cubature::cuhre     |   642µs |   666µs |  1495.135 |   39.84KB |  30.513 |    49 |    1 |    32.77ms |

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
| Non-vectorized Hcubature       | 1.35ms | 1.39ms | 715.864 |    85.2KB | 37.677 |    19 |    1 |     26.5ms |
| Vectorized Hcubature           | 1.78ms | 1.83ms | 547.057 |   147.2KB |  0.000 |    20 |    0 |     36.6ms |
| Non-vectorized Pcubature       | 2.02ms | 2.03ms | 484.986 |        0B | 53.887 |    18 |    2 |     37.1ms |
| Vectorized Pcubature           | 2.56ms |  2.6ms | 381.754 |    68.9KB | 20.092 |    19 |    1 |     49.8ms |
| Non-vectorized cubature::cuhre | 3.22ms | 3.27ms | 304.090 |        0B | 33.788 |    18 |    2 |     59.2ms |
| Vectorized cubature::cuhre     | 3.71ms | 3.77ms | 264.152 |   125.6KB | 13.903 |    19 |    1 |     71.9ms |

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
| Non-vectorized Hcubature       | 3.58ms | 3.63ms | 273.069 |     133KB | 48.189 |    17 |    3 |     62.3ms |
| Vectorized Hcubature           | 4.58ms | 4.64ms | 214.781 |     249KB | 23.865 |    18 |    2 |     83.8ms |
| Non-vectorized Pcubature       |  2.5ms | 2.53ms | 391.725 |        0B | 43.525 |    18 |    2 |       46ms |
| Vectorized Pcubature           | 3.23ms |  3.3ms | 301.751 |      69KB | 15.882 |    19 |    1 |       63ms |
| Non-vectorized cubature::cuhre | 6.94ms | 7.01ms | 142.227 |        0B | 35.557 |    16 |    4 |    112.5ms |
| Vectorized cubature::cuhre     | 8.33ms |  8.4ms | 118.700 |     224KB | 29.675 |    16 |    4 |    134.8ms |

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

| expression                     |    min |  median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|-------:|--------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 1.59ms |  1.62ms | 616.056 |    64.6KB | 32.424 |    19 |    1 |     30.8ms |
| Vectorized Hcubature           | 2.11ms |  2.19ms | 458.246 |     156KB | 24.118 |    19 |    1 |     41.5ms |
| Non-vectorized Pcubature       | 8.15ms |  8.33ms | 120.320 |        0B | 40.107 |    15 |    5 |    124.7ms |
| Vectorized Pcubature           | 10.3ms | 10.43ms |  95.765 |   386.2KB | 23.941 |    16 |    4 |    167.1ms |
| Non-vectorized cubature::cuhre | 4.32ms |  4.41ms | 226.213 |        0B | 39.920 |    17 |    3 |     75.2ms |
| Vectorized cubature::cuhre     | 4.69ms |  4.82ms | 201.774 |   225.8KB | 22.419 |    18 |    2 |     89.2ms |

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
| Non-vectorized Hcubature       |  3.32ms |  3.39ms | 278.260 |   32.89KB | 27.826 |    20 |    2 |     71.9ms |
| Vectorized Hcubature           |  3.96ms |  4.11ms | 233.202 |  205.61KB | 23.320 |    20 |    2 |     85.8ms |
| Non-vectorized Pcubature       |  8.25ms |  8.41ms | 109.912 |        0B | 32.974 |    20 |    6 |      182ms |
| Vectorized Pcubature           |  9.47ms |  9.66ms |  97.057 |  386.24KB | 24.264 |    20 |    5 |    206.1ms |
| Non-vectorized cubature::cuhre | 42.52ms |  44.2ms |  22.717 |        0B | 26.124 |    20 |   23 |    880.4ms |
| Vectorized cubature::cuhre     | 42.82ms | 43.41ms |  22.817 |    2.04MB | 25.099 |    20 |   22 |    876.5ms |

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
| Non-vectorized Hcubature       | 128.6µs | 133.1µs |  7343.472 |    53.7KB |  0.000 |   100 |    0 |    13.62ms |
| Vectorized Hcubature           | 210.7µs | 221.5µs |  4361.072 |    68.5KB | 44.051 |    99 |    1 |     22.7ms |
| Non-vectorized Pcubature       |  54.7µs |  57.7µs | 16599.833 |        0B |  0.000 |   100 |    0 |     6.02ms |
| Vectorized Pcubature           |   157µs | 162.9µs |  5918.349 |        0B |  0.000 |   100 |    0 |     16.9ms |
| Non-vectorized cubature::cuhre | 222.6µs | 230.2µs |  4254.330 |        0B | 42.973 |    99 |    1 |    23.27ms |
| Vectorized cubature::cuhre     | 498.5µs | 527.6µs |  1880.729 |        0B | 18.997 |    99 |    1 |    52.64ms |

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
| Non-vectorized Hcubature       |    4.5ms |   4.58ms |  217.292 |    78.4KB | 35.373 |    86 |   14 |    395.8ms |
| Vectorized Hcubature           |   5.32ms |    5.4ms |  184.913 |   304.7KB | 27.631 |    87 |   13 |    470.5ms |
| Non-vectorized Pcubature       | 409.23µs | 426.06µs | 2334.308 |        0B | 23.579 |    99 |    1 |     42.4ms |
| Vectorized Pcubature           | 599.31µs | 634.03µs | 1566.569 |    18.3KB | 31.971 |    98 |    2 |     62.6ms |
| Non-vectorized cubature::cuhre |   1.21ms |   1.24ms |  799.006 |        0B | 24.712 |    97 |    3 |    121.4ms |
| Vectorized cubature::cuhre     |   1.32ms |   1.35ms |  734.291 |    60.1KB | 22.710 |    97 |    3 |    132.1ms |

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
