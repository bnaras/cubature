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
| Non-vectorized Hcubature       | 271.15µs | 282.5µs | 3500.916 |    59.3KB | 35.363 |    99 |    1 |     28.3ms |
| Vectorized Hcubature           | 385.87µs | 409.5µs | 2407.857 |    54.5KB | 24.322 |    99 |    1 |     41.1ms |
| Non-vectorized Pcubature       | 853.39µs | 877.5µs | 1121.619 |    31.2KB | 34.689 |    97 |    3 |     86.5ms |
| Vectorized Pcubature           |   1.17ms |   1.2ms |  825.996 |    58.4KB | 16.857 |    98 |    2 |    118.6ms |
| Non-vectorized cubature::cuhre | 582.81µs | 599.6µs | 1652.723 |    40.5KB | 16.694 |    99 |    1 |     59.9ms |
| Vectorized cubature::cuhre     | 591.27µs | 626.2µs | 1594.027 |    39.8KB | 32.531 |    98 |    2 |     61.5ms |

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
| Non-vectorized Hcubature       | 815.46ms | 823.41ms |   1.192 |  139.68KB | 18.477 |    10 |  155 |      8.39s |
| Vectorized Hcubature           |   1.79ms |   2.53ms | 427.225 |    1.89MB |  0.000 |    10 |    0 |    23.41ms |
| Non-vectorized Pcubature       | 349.74ms | 352.99ms |   2.790 |        0B | 18.137 |    10 |   65 |      3.58s |
| Vectorized Pcubature           |   1.26ms |   1.28ms | 662.610 |  810.26KB | 66.261 |    10 |    1 |    15.09ms |
| Non-vectorized cubature::cuhre | 335.21ms | 338.69ms |   2.939 |        0B | 18.514 |    10 |   63 |       3.4s |
| Vectorized cubature::cuhre     |   3.23ms |   3.29ms | 302.477 |  898.41KB |  0.000 |    10 |    0 |    33.06ms |

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
| g1()       | 755µs | 1.37ms | 732.882 |  355.13KB |      0 |    20 |    0 |     27.3ms |
| g2()       | 753µs | 1.37ms | 758.122 |    2.49KB |      0 |    20 |    0 |     26.4ms |
| g3()       | 749µs | 1.37ms | 733.728 |    2.49KB |      0 |    20 |    0 |     27.3ms |

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
| Non-vectorized Hcubature       |  36.8µs |  38.9µs | 24897.671 |    7.66KB | 24.923 |   999 |    1 |     40.1ms |
| Vectorized Hcubature           |    64µs |  67.5µs | 14392.020 |    1.16KB | 28.842 |   998 |    2 |     69.3ms |
| Non-vectorized Pcubature       |    49µs |    51µs | 19015.257 |        0B | 19.034 |   999 |    1 |     52.5ms |
| Vectorized Pcubature           | 113.1µs | 118.3µs |  8231.283 |   18.68KB | 24.768 |   997 |    3 |    121.1ms |
| Non-vectorized cubature::cuhre | 348.6µs | 360.4µs |  2760.045 |        0B | 25.066 |   991 |    9 |    359.1ms |
| Vectorized cubature::cuhre     | 363.5µs | 380.7µs |  2608.107 |   16.38KB | 26.345 |   990 |   10 |    379.6ms |

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
| Non-vectorized Hcubature       |   2.92ms |   2.94ms |   338.021 |    67.3KB |  0.000 |    10 |    0 |    29.58ms |
| Vectorized Hcubature           |      5ms |   5.08ms |   196.577 |   290.5KB | 49.144 |     8 |    2 |     40.7ms |
| Non-vectorized Pcubature       |  74.94µs |  75.69µs | 12142.883 |        0B |  0.000 |    10 |    0 |   823.53µs |
| Vectorized Pcubature           | 152.34µs | 156.76µs |  6169.471 |     4.1KB |  0.000 |    10 |    0 |     1.62ms |
| Non-vectorized cubature::cuhre |  13.47ms |  13.62ms |    73.389 |        0B | 31.452 |     7 |    3 |    95.38ms |
| Vectorized cubature::cuhre     |   20.4ms |  20.57ms |    47.606 |   971.5KB | 47.606 |     5 |    5 |   105.03ms |

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
| Non-vectorized Hcubature       |    45ms |  45.7ms |  18.820 |    17.7KB | 37.640 |    10 |   20 |   531.35ms |
| Vectorized Hcubature           |  48.3ms |    49ms |  20.193 |  1011.8KB | 26.251 |    10 |   13 |   495.23ms |
| Non-vectorized cubature::cuhre | 821.7ms | 831.1ms |   1.198 |        0B | 28.158 |    10 |  235 |      8.35s |
| Vectorized cubature::cuhre     | 901.8ms | 910.5ms |   1.087 |    21.2MB | 24.682 |    10 |  227 |       9.2s |

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
| Non-vectorized Hcubature       |  57.4µs |  59.4µs | 16111.269 |    6.75KB |   0.000 |    50 |    0 |      3.1ms |
| Vectorized Hcubature           |  91.9µs |  95.1µs | 10006.489 |    2.91KB |   0.000 |    50 |    0 |        5ms |
| Non-vectorized Pcubature       |  52.9µs |  54.4µs | 17400.128 |        0B | 355.105 |    49 |    1 |     2.82ms |
| Vectorized Pcubature           |  81.5µs |  83.5µs | 11536.331 |   19.12KB |   0.000 |    50 |    0 |     4.33ms |
| Non-vectorized cubature::cuhre | 658.4µs | 681.4µs |  1467.779 |        0B |  29.955 |    49 |    1 |    33.38ms |
| Vectorized cubature::cuhre     | 658.3µs | 688.2µs |  1444.840 |   39.84KB |   0.000 |    50 |    0 |    34.61ms |

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
| Non-vectorized Hcubature       | 1.36ms |  1.4ms | 712.594 |    85.2KB | 37.505 |    19 |    1 |     26.7ms |
| Vectorized Hcubature           | 1.79ms | 1.83ms | 540.368 |   147.2KB | 28.440 |    19 |    1 |     35.2ms |
| Non-vectorized Pcubature       | 2.04ms | 2.08ms | 478.997 |        0B | 25.210 |    19 |    1 |     39.7ms |
| Vectorized Pcubature           |  2.6ms | 2.68ms | 372.457 |    68.9KB | 41.384 |    18 |    2 |     48.3ms |
| Non-vectorized cubature::cuhre | 3.27ms | 3.32ms | 300.503 |        0B | 33.389 |    18 |    2 |     59.9ms |
| Vectorized cubature::cuhre     | 3.75ms | 3.84ms | 259.302 |   125.6KB | 13.647 |    19 |    1 |     73.3ms |

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
| Non-vectorized Hcubature       | 3.61ms |  3.7ms | 269.823 |     133KB | 47.616 |    17 |    3 |       63ms |
| Vectorized Hcubature           | 4.64ms | 4.75ms | 209.857 |     249KB | 23.317 |    18 |    2 |     85.8ms |
| Non-vectorized Pcubature       | 2.51ms | 2.53ms | 393.524 |        0B | 43.725 |    18 |    2 |     45.7ms |
| Vectorized Pcubature           | 3.31ms | 3.38ms | 296.237 |      69KB | 15.591 |    19 |    1 |     64.1ms |
| Non-vectorized cubature::cuhre | 6.94ms | 7.07ms | 141.129 |        0B | 47.043 |    15 |    5 |    106.3ms |
| Vectorized cubature::cuhre     | 8.29ms | 8.34ms | 119.569 |     224KB | 21.100 |    17 |    3 |    142.2ms |

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
| Non-vectorized Hcubature       |  1.57ms |   1.6ms | 625.195 |    64.6KB | 32.905 |    19 |    1 |     30.4ms |
| Vectorized Hcubature           |  2.14ms |  2.18ms | 458.268 |     156KB | 24.119 |    19 |    1 |     41.5ms |
| Non-vectorized Pcubature       |  8.14ms |  8.24ms | 121.362 |        0B | 40.454 |    15 |    5 |    123.6ms |
| Vectorized Pcubature           | 10.36ms | 10.51ms |  95.264 |   386.2KB | 31.755 |    15 |    5 |    157.5ms |
| Non-vectorized cubature::cuhre |  4.43ms |   4.5ms | 221.921 |        0B | 24.658 |    18 |    2 |     81.1ms |
| Vectorized cubature::cuhre     |  4.81ms |  4.91ms | 203.462 |   225.8KB | 22.607 |    18 |    2 |     88.5ms |

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
| Non-vectorized Hcubature       |  3.34ms |  3.43ms | 276.057 |   32.89KB | 27.606 |    20 |    2 |     72.4ms |
| Vectorized Hcubature           |  4.12ms |  4.23ms | 224.762 |  205.61KB | 22.476 |    20 |    2 |       89ms |
| Non-vectorized Pcubature       |  8.32ms |  8.57ms | 111.060 |        0B | 27.765 |    20 |    5 |    180.1ms |
| Vectorized Pcubature           |  9.62ms |  9.89ms |  95.935 |  386.24KB | 23.984 |    20 |    5 |    208.5ms |
| Non-vectorized cubature::cuhre | 43.63ms | 44.13ms |  22.428 |        0B | 24.670 |    20 |   22 |    891.8ms |
| Vectorized cubature::cuhre     | 43.45ms | 44.21ms |  22.497 |    2.04MB | 23.622 |    20 |   21 |      889ms |

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

| expression                     |   min | median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|------:|-------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 126µs |  135µs |  7085.618 |    53.7KB | 71.572 |    99 |    1 |    13.97ms |
| Vectorized Hcubature           | 210µs |  220µs |  4421.691 |    68.5KB |  0.000 |   100 |    0 |    22.62ms |
| Non-vectorized Pcubature       |  52µs |   54µs | 17799.208 |        0B |  0.000 |   100 |    0 |     5.62ms |
| Vectorized Pcubature           | 157µs |  164µs |  5752.952 |        0B | 58.111 |    99 |    1 |    17.21ms |
| Non-vectorized cubature::cuhre | 226µs |  233µs |  4241.178 |        0B |  0.000 |   100 |    0 |    23.58ms |
| Vectorized cubature::cuhre     | 507µs |  527µs |  1869.567 |        0B | 18.885 |    99 |    1 |    52.95ms |

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
| Non-vectorized Hcubature       |   4.54ms |   4.62ms |  215.203 |    78.4KB | 35.033 |    86 |   14 |    399.6ms |
| Vectorized Hcubature           |   5.41ms |   5.55ms |  179.898 |   304.7KB | 26.881 |    87 |   13 |    483.6ms |
| Non-vectorized Pcubature       | 406.38µs | 420.48µs | 2370.203 |        0B | 23.941 |    99 |    1 |     41.8ms |
| Vectorized Pcubature           | 600.28µs |  627.3µs | 1574.371 |    18.3KB | 32.130 |    98 |    2 |     62.2ms |
| Non-vectorized cubature::cuhre |   1.24ms |   1.26ms |  793.253 |        0B | 24.534 |    97 |    3 |    122.3ms |
| Vectorized cubature::cuhre     |   1.33ms |   1.37ms |  721.827 |    60.1KB | 22.325 |    97 |    3 |    134.4ms |

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
    ## [25] digest_0.6.39     R6_2.6.1          pillar_1.11.1     magrittr_2.0.4   
    ## [29] bslib_0.9.0       tools_4.5.2       pkgdown_2.2.0     cachem_1.1.0     
    ## [33] desc_1.4.3
