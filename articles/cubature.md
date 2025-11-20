# Cubature Vectorization Results

## Introduction

This R `cubature` package exposes both the `hcubature` and `pcubature`
routines of the underlying C `cubature` library, including the
vectorized interfaces.

Per the documentation, use of `pcubature` is advisable only for smooth
integrands in dimensions up to three at most. In fact, the `pcubature`
routines perform significantly worse than the vectorized `hcubature` in
inappropriate cases. So when in doubt, you are better off using
`hcubature`.

Version 2.0 of this package integrates the `Cuba` library as well, once
again providing vectorized interfaces.

The main point of this note is to examine the difference vectorization
makes. My recommendations are below in the summary section.

## A Timing Harness

Our harness will provide timing results for `hcubature`, `pcubature`
(where appropriate) and Cuba `cuhre` calls. We begin by creating a
harness for these calls.

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

## Example 1.

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

| expression                     |      min |  median |  itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|--------:|---------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 273.48µs |   285µs | 3490.204 |    59.3KB | 35.255 |    99 |    1 |     28.4ms |
| Vectorized Hcubature           | 390.58µs | 410.9µs | 2401.619 |    54.5KB | 24.259 |    99 |    1 |     41.2ms |
| Non-vectorized Pcubature       | 860.38µs | 884.7µs | 1108.807 |    31.2KB | 34.293 |    97 |    3 |     87.5ms |
| Vectorized Pcubature           |   1.17ms |   1.2ms |  826.518 |    58.4KB | 16.868 |    98 |    2 |    118.6ms |
| Non-vectorized cubature::cuhre | 579.96µs | 602.2µs | 1646.980 |    40.5KB | 16.636 |    99 |    1 |     60.1ms |
| Vectorized cubature::cuhre     | 592.14µs | 627.2µs | 1593.192 |    39.8KB | 32.514 |    98 |    2 |     61.5ms |

## Multivariate Normal

Using `cubature`, we evaluate $$\int_{R}\phi(x)dx$$ where $\phi(x)$ is
the three-dimensional multivariate normal density with mean 0, and
variance $$\Sigma = \left( \begin{array}{rrr}
1 & \frac{3}{5} & \frac{1}{3} \\
\frac{3}{5} & 1 & \frac{11}{15} \\
\frac{1}{3} & \frac{11}{15} & 1
\end{array} \right)$$ and $R$ is
$\left\lbrack - \frac{1}{2},1 \right\rbrack \times \left\lbrack - \frac{1}{2},4 \right\rbrack \times \left\lbrack - \frac{1}{2},2 \right\rbrack.$

We construct a scalar function (`my_dmvnorm`) and a vector analog
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
| Non-vectorized Hcubature       | 803.24ms | 819.11ms |   1.209 |  139.68KB | 18.743 |    10 |  155 |      8.27s |
| Vectorized Hcubature           |    1.8ms |   2.41ms | 441.466 |    1.89MB |  0.000 |    10 |    0 |    22.65ms |
| Non-vectorized Pcubature       | 346.98ms | 350.18ms |   2.813 |        0B | 18.286 |    10 |   65 |      3.56s |
| Vectorized Pcubature           |   1.24ms |   1.26ms | 672.872 |  810.26KB | 67.287 |    10 |    1 |    14.86ms |
| Non-vectorized cubature::cuhre |  336.2ms | 341.02ms |   2.939 |        0B | 18.514 |    10 |   63 |       3.4s |
| Vectorized cubature::cuhre     |    3.2ms |   3.25ms | 306.143 |  898.41KB |  0.000 |    10 |    0 |    32.66ms |

The effect of vectorization is huge. So it makes sense for users to
vectorize the integrands as much as possible for efficiency.

Furthermore, for this particular example, we expect
[`mvtnorm::pmvnorm`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html) to do
pretty well since it is specialized for the multivariate normal. The
vectorized versions of `hcubature` and `pcubature` seem competitive and
in some cases better, if you compare the table above to the one below.

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
| g1()       | 771µs | 1.37ms | 737.686 |  355.13KB |      0 |    20 |    0 |     27.1ms |
| g2()       | 754µs | 1.37ms | 757.994 |    2.49KB |      0 |    20 |    0 |     26.4ms |
| g3()       | 762µs | 1.37ms | 741.815 |    2.49KB |      0 |    20 |    0 |       27ms |

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
| Non-vectorized Hcubature       |  37.1µs |  39.6µs | 24537.208 |    7.66KB | 24.562 |   999 |    1 |     40.7ms |
| Vectorized Hcubature           |  65.1µs |  68.4µs | 14096.474 |    1.16KB | 28.249 |   998 |    2 |     70.8ms |
| Non-vectorized Pcubature       |  48.8µs |  50.8µs | 19151.815 |        0B | 19.171 |   999 |    1 |     52.2ms |
| Vectorized Pcubature           | 113.9µs | 120.6µs |  8126.376 |   18.68KB | 24.452 |   997 |    3 |    122.7ms |
| Non-vectorized cubature::cuhre | 331.9µs | 342.4µs |  2907.356 |        0B | 26.404 |   991 |    9 |    340.9ms |
| Vectorized cubature::cuhre     | 360.2µs | 377.4µs |  2637.858 |   16.38KB | 26.645 |   990 |   10 |    375.3ms |

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
| Non-vectorized Hcubature       |   2.85ms |   2.91ms |   343.528 |    67.3KB |  0.000 |    10 |    0 |    29.11ms |
| Vectorized Hcubature           |   4.88ms |   4.99ms |   200.485 |   290.5KB | 22.276 |     9 |    1 |    44.89ms |
| Non-vectorized Pcubature       |  73.12µs |  74.69µs | 12247.688 |        0B |  0.000 |    10 |    0 |   816.48µs |
| Vectorized Pcubature           | 152.36µs | 155.59µs |  6200.883 |     4.1KB |  0.000 |    10 |    0 |     1.61ms |
| Non-vectorized cubature::cuhre |  13.43ms |   13.6ms |    73.776 |        0B | 31.618 |     7 |    3 |    94.88ms |
| Vectorized cubature::cuhre     |  20.14ms |   20.2ms |    49.389 |   971.5KB | 49.389 |     5 |    5 |   101.24ms |

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
| Non-vectorized Hcubature       |  44.2ms |  45.4ms |  18.880 |    17.7KB | 37.760 |    10 |   20 |   529.66ms |
| Vectorized Hcubature           |  47.7ms |  48.6ms |  20.469 |  1011.8KB | 26.609 |    10 |   13 |   488.55ms |
| Non-vectorized cubature::cuhre | 794.6ms | 798.8ms |   1.240 |        0B | 29.129 |    10 |  235 |      8.07s |
| Vectorized cubature::cuhre     | 884.2ms | 900.5ms |   1.100 |    21.2MB | 24.980 |    10 |  227 |      9.09s |

## A Simple Polynomial (product of coordinates)

``` r
testFn3 <- function(x) prod(2 * x)
testFn3_v <- function(x) matrix(apply(x, 2, function(z) prod(2 * z)), ncol = ncol(x))

d <- harness(f = testFn3, fv = testFn3_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 50)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median |   itr/sec | mem_alloc |  gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|----------:|----------:|--------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |    57µs |  59.4µs | 16165.859 |    6.75KB |   0.000 |    50 |    0 |     3.09ms |
| Vectorized Hcubature           |  89.8µs |  94.8µs | 10118.425 |    2.91KB |   0.000 |    50 |    0 |     4.94ms |
| Non-vectorized Pcubature       |  51.5µs |  54.2µs | 17801.281 |        0B |   0.000 |    50 |    0 |     2.81ms |
| Vectorized Pcubature           |  80.5µs |  83.3µs | 11102.429 |   19.12KB | 226.580 |    49 |    1 |     4.41ms |
| Non-vectorized cubature::cuhre | 634.6µs |   657µs |  1514.700 |        0B |  30.912 |    49 |    1 |    32.35ms |
| Vectorized cubature::cuhre     | 640.4µs | 669.4µs |  1466.556 |   39.84KB |   0.000 |    50 |    0 |    34.09ms |

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
| Non-vectorized Hcubature       | 1.34ms | 1.38ms | 725.133 |    85.2KB | 38.165 |    19 |    1 |     26.2ms |
| Vectorized Hcubature           | 1.78ms |  1.8ms | 546.771 |   147.2KB | 28.777 |    19 |    1 |     34.7ms |
| Non-vectorized Pcubature       | 1.99ms | 2.03ms | 490.126 |        0B | 25.796 |    19 |    1 |     38.8ms |
| Vectorized Pcubature           | 2.57ms | 2.64ms | 378.722 |    68.9KB | 42.080 |    18 |    2 |     47.5ms |
| Non-vectorized cubature::cuhre |  3.2ms | 3.26ms | 305.959 |        0B | 33.995 |    18 |    2 |     58.8ms |
| Vectorized cubature::cuhre     | 3.71ms | 3.79ms | 263.607 |   125.6KB | 13.874 |    19 |    1 |     72.1ms |

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
| Non-vectorized Hcubature       | 3.59ms | 3.65ms | 273.151 |     133KB | 48.203 |    17 |    3 |     62.2ms |
| Vectorized Hcubature           | 4.67ms | 4.73ms | 210.235 |     249KB | 23.359 |    18 |    2 |     85.6ms |
| Non-vectorized Pcubature       |  2.5ms | 2.53ms | 394.095 |        0B | 43.788 |    18 |    2 |     45.7ms |
| Vectorized Pcubature           | 3.27ms | 3.33ms | 300.829 |      69KB | 15.833 |    19 |    1 |     63.2ms |
| Non-vectorized cubature::cuhre | 6.94ms | 7.05ms | 141.989 |        0B | 47.330 |    15 |    5 |    105.6ms |
| Vectorized cubature::cuhre     | 8.22ms | 8.31ms | 120.213 |     224KB | 21.214 |    17 |    3 |    141.4ms |

## Tsuda’s Example

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
| Non-vectorized Hcubature       |  1.58ms |  1.59ms | 625.622 |    64.6KB |  0.000 |    20 |    0 |       32ms |
| Vectorized Hcubature           |   2.1ms |  2.13ms | 468.749 |     156KB | 24.671 |    19 |    1 |     40.5ms |
| Non-vectorized Pcubature       |  8.06ms |  8.16ms | 122.093 |        0B | 40.698 |    15 |    5 |    122.9ms |
| Vectorized Pcubature           | 10.25ms | 10.41ms |  96.101 |   386.2KB | 32.034 |    15 |    5 |    156.1ms |
| Non-vectorized cubature::cuhre |  4.27ms |  4.35ms | 230.083 |        0B | 25.565 |    18 |    2 |     78.2ms |
| Vectorized cubature::cuhre     |  4.77ms |  4.87ms | 204.768 |   225.8KB | 22.752 |    18 |    2 |     87.9ms |

## Morokoff & Calflish Example

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
| Non-vectorized Hcubature       |   3.3ms |  3.41ms | 271.890 |   32.89KB | 27.189 |    20 |    2 |     73.6ms |
| Vectorized Hcubature           |  4.04ms |  4.16ms | 226.515 |  205.61KB | 22.652 |    20 |    2 |     88.3ms |
| Non-vectorized Pcubature       |  8.27ms |  8.57ms | 110.751 |        0B | 27.688 |    20 |    5 |    180.6ms |
| Vectorized Pcubature           |   9.5ms |  9.72ms |  97.532 |  386.24KB | 24.383 |    20 |    5 |    205.1ms |
| Non-vectorized cubature::cuhre | 43.03ms | 43.81ms |  22.604 |        0B | 24.865 |    20 |   22 |    884.8ms |
| Vectorized cubature::cuhre     | 42.86ms | 43.69ms |  22.767 |    2.04MB | 23.905 |    20 |   21 |    878.5ms |

## Wang-Landau Sampling 1d, 2d Examples

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
| Non-vectorized Hcubature       | 126.9µs |   137µs |  7078.559 |    53.7KB | 71.501 |    99 |    1 |       14ms |
| Vectorized Hcubature           | 210.3µs | 220.8µs |  4416.320 |    68.5KB |  0.000 |   100 |    0 |     22.6ms |
| Non-vectorized Pcubature       |  51.6µs |  54.2µs | 17867.145 |        0B |  0.000 |   100 |    0 |      5.6ms |
| Vectorized Pcubature           | 155.2µs | 165.5µs |  5926.995 |        0B | 59.869 |    99 |    1 |     16.7ms |
| Non-vectorized cubature::cuhre | 220.7µs | 228.5µs |  4309.506 |        0B |  0.000 |   100 |    0 |     23.2ms |
| Vectorized cubature::cuhre     | 493.7µs | 516.8µs |  1923.373 |        0B | 19.428 |    99 |    1 |     51.5ms |

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
| Non-vectorized Hcubature       |   4.55ms |   4.62ms |  215.095 |    78.4KB | 35.015 |    86 |   14 |    399.8ms |
| Vectorized Hcubature           |   5.39ms |    5.5ms |  181.368 |   304.7KB | 27.101 |    87 |   13 |    479.7ms |
| Non-vectorized Pcubature       | 409.65µs | 426.18µs | 2335.077 |        0B | 23.587 |    99 |    1 |     42.4ms |
| Vectorized Pcubature           | 600.66µs | 628.89µs | 1577.618 |    18.3KB | 32.196 |    98 |    2 |     62.1ms |
| Non-vectorized cubature::cuhre |   1.23ms |   1.25ms |  788.323 |        0B | 24.381 |    97 |    3 |      123ms |
| Vectorized cubature::cuhre     |   1.32ms |   1.36ms |  726.245 |    60.1KB | 22.461 |    97 |    3 |    133.6ms |

## Implementation Notes

About the only real modification we have made to the underlying
`cubature` library is that we use `M = 16` rather than the default
`M = 19` suggested by the original author for `pcubature`. This allows
us to comply with CRAN package size limits and seems to work reasonably
well for the above tests. Future versions will allow for such
customization on demand.

The changes made to the `Cuba` library are managed in a [Github repo
branch](https://github.com/bnaras/Cuba/tree/R_pkg): each time a new
release is made, we update the main branch, and keep all changes for
Unix platforms in a branch named `R_pkg` against the current main
branch. Customization for windows is done in the package itself using
the `Makevars.win` script.

## Summary

The recommendations are:

1.  *Vectorize* your function. The time spent in so doing pays back
    enormously. This is easy to do and the examples above show how.

2.  Vectorized `hcubature` seems to be a good starting point.

3.  For smooth integrands in low dimensions ($\leq 3$), `pcubature`
    might be worth trying out. Experiment before using in a production
    package.

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
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
    ## [1] mvtnorm_1.3-3    cubature_2.1.4-3 bench_1.1.4     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.6.5       cli_3.6.5         knitr_1.50        rlang_1.1.6      
    ##  [5] xfun_0.54         textshaping_1.0.4 jsonlite_2.0.0    glue_1.8.0       
    ##  [9] htmltools_0.5.8.1 ragg_1.5.0        sass_0.4.10       rmarkdown_2.30   
    ## [13] tibble_3.3.0      evaluate_1.0.5    jquerylib_0.1.4   fastmap_1.2.0    
    ## [17] profmem_0.7.0     yaml_2.3.10       lifecycle_1.0.4   compiler_4.5.2   
    ## [21] fs_1.6.6          pkgconfig_2.0.3   Rcpp_1.1.0        systemfonts_1.3.1
    ## [25] digest_0.6.38     R6_2.6.1          pillar_1.11.1     magrittr_2.0.4   
    ## [29] bslib_0.9.0       tools_4.5.2       pkgdown_2.2.0     cachem_1.1.0     
    ## [33] desc_1.4.3
