# Version 2.0 Notes

## What’s new

Version 2.0 integrates two well-known cubature libraries in one place:

- The [cubature C library](https://github.com/stevengj/cubature) of
  Steven G. Johnson.
- The [Cuba C library](https://feynarts.de/cuba/) of Thomas Hahn.

It also provides a single function `cubintegrate` that allows one to
call all methods in a uniform fashion, as I explain below.

**N.B.** One has to be aware that there are cases where one library will
integrate a function while the other won’t, and in some cases, provide
somewhat different answers. That still makes sense and depends on the
underlying methodology used.

## Unified Interface

Following a suggestion by Simen Guare, we now have a function
`cubintegrate` that can be used to try out various integration methods
easily. Some examples.

``` r
library(cubature)
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
```

First we try the scalar invocation with `hcubature`.

``` r
cubintegrate(f = my_dmvnorm, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "pcubature",
             mean = rep(0, m), sigma = sigma, logdet = logdet)
```

    ## $integral
    ## [1] 0.3341125
    ## 
    ## $error
    ## [1] 2.21207e-06
    ## 
    ## $neval
    ## [1] 4913
    ## 
    ## $returnCode
    ## [1] 0

We can compare that with Cuba’s `cuhre`.

``` r
cubintegrate(f = my_dmvnorm, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "cuhre",
             mean = rep(0, m), sigma = sigma, logdet = logdet)
```

    ## $integral
    ## [1] 0.3341125
    ## 
    ## $error
    ## [1] 2.226266e-06
    ## 
    ## $nregions
    ## [1] 19
    ## 
    ## $neval
    ## [1] 4699
    ## 
    ## $prob
    ## [1] 0
    ## 
    ## $returnCode
    ## [1] 0

The Cuba routine can take various further arguments; see for example,
the help on `cuhre`. Such arguments can be directly passed to
`cubintegrate`.

``` r
cubintegrate(f = my_dmvnorm, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "cuhre",
             mean = rep(0, m), sigma = sigma, logdet = logdet,
             flags = list(verbose = 2))
```

    ## Cuhre input parameters:
    ##   ndim 3
    ##   ncomp 1
    ##   nvec 1
    ##   epsrel 1e-05
    ##   epsabs 1e-12
    ##   flags 6
    ##   mineval 0
    ##   maxeval 1000000
    ##   key 0
    ##   statefile "(null)"
    ## 
    ## Iteration 1:  127 integrand evaluations so far
    ## [1] 0.330292 +- 0.16199      chisq 0 (0 df)
    ## 
    ## Iteration 2:  381 integrand evaluations so far
    ## [1] 0.334551 +- 0.0726708    chisq 0.000575398 (1 df)
    ## 
    ## Iteration 3:  635 integrand evaluations so far
    ## [1] 0.334084 +- 0.0126583    chisq 0.000588724 (2 df)
    ## 
    ## Iteration 4:  889 integrand evaluations so far
    ## [1] 0.334093 +- 0.00467064   chisq 0.000590497 (3 df)
    ## 
    ## Iteration 5:  1143 integrand evaluations so far
    ## [1] 0.33411 +- 0.00213905    chisq 0.00060536 (4 df)
    ## 
    ## Iteration 6:  1397 integrand evaluations so far
    ## [1] 0.334113 +- 0.000658593      chisq 0.000616501 (5 df)
    ## 
    ## Iteration 7:  1651 integrand evaluations so far
    ## [1] 0.334113 +- 0.000290809      chisq 0.000617561 (6 df)
    ## 
    ## Iteration 8:  1905 integrand evaluations so far
    ## [1] 0.334113 +- 4.6551e-05   chisq 0.000617624 (7 df)
    ## 
    ## Iteration 9:  2159 integrand evaluations so far
    ## [1] 0.334113 +- 3.69574e-05      chisq 0.000617628 (8 df)
    ## 
    ## Iteration 10:  2413 integrand evaluations so far
    ## [1] 0.334113 +- 6.66791e-05      chisq 0.000623039 (9 df)
    ## 
    ## Iteration 11:  2667 integrand evaluations so far
    ## [1] 0.334113 +- 2.65953e-05      chisq 0.000637044 (10 df)
    ## 
    ## Iteration 12:  2921 integrand evaluations so far
    ## [1] 0.334113 +- 3.88133e-05      chisq 0.000639547 (11 df)
    ## 
    ## Iteration 13:  3175 integrand evaluations so far
    ## [1] 0.334113 +- 2.10511e-05      chisq 0.000643649 (12 df)
    ## 
    ## Iteration 14:  3429 integrand evaluations so far
    ## [1] 0.334113 +- 1.31307e-05      chisq 0.000650985 (13 df)
    ## 
    ## Iteration 15:  3683 integrand evaluations so far
    ## [1] 0.334113 +- 8.69938e-06      chisq 0.000654884 (14 df)
    ## 
    ## Iteration 16:  3937 integrand evaluations so far
    ## [1] 0.334113 +- 1.88417e-05      chisq 0.000655659 (15 df)
    ## 
    ## Iteration 17:  4191 integrand evaluations so far
    ## [1] 0.334113 +- 1.59906e-05      chisq 0.000656612 (16 df)
    ## 
    ## Iteration 18:  4445 integrand evaluations so far
    ## [1] 0.334113 +- 4.85581e-06      chisq 0.000660539 (17 df)
    ## 
    ## Iteration 19:  4699 integrand evaluations so far
    ## [1] 0.334113 +- 2.22627e-06      chisq 0.000668228 (18 df)

    ## $integral
    ## [1] 0.3341125
    ## 
    ## $error
    ## [1] 2.226266e-06
    ## 
    ## $nregions
    ## [1] 19
    ## 
    ## $neval
    ## [1] 4699
    ## 
    ## $prob
    ## [1] 0
    ## 
    ## $returnCode
    ## [1] 0

As there are many such method-specific arguments, you may find the
function [`default_args()`](../reference/default_args.md) useful.

``` r
str(default_args())
```

    ## List of 6
    ##  $ hcubature:List of 1
    ##   ..$ norm: chr [1:5] "INDIVIDUAL" "PAIRED" "L2" "L1" ...
    ##  $ pcubature:List of 1
    ##   ..$ norm: chr [1:5] "INDIVIDUAL" "PAIRED" "L2" "L1" ...
    ##  $ cuhre    :List of 4
    ##   ..$ minEval  : int 0
    ##   ..$ stateFile: NULL
    ##   ..$ flags    :List of 6
    ##   .. ..$ verbose   : int 0
    ##   .. ..$ final     : int 1
    ##   .. ..$ smooth    : int 0
    ##   .. ..$ keep_state: int 0
    ##   .. ..$ load_state: int 0
    ##   .. ..$ level     : int 0
    ##   ..$ key      : int 0
    ##  $ divonne  :List of 14
    ##   ..$ minEval     : int 0
    ##   ..$ stateFile   : NULL
    ##   ..$ flags       :List of 6
    ##   .. ..$ verbose   : int 0
    ##   .. ..$ final     : int 1
    ##   .. ..$ smooth    : int 0
    ##   .. ..$ keep_state: int 0
    ##   .. ..$ load_state: int 0
    ##   .. ..$ level     : int 0
    ##   ..$ rngSeed     : int 0
    ##   ..$ key1        : int 47
    ##   ..$ key2        : int 1
    ##   ..$ key3        : int 1
    ##   ..$ maxPass     : int 5
    ##   ..$ border      : num 0
    ##   ..$ maxChisq    : num 10
    ##   ..$ minDeviation: num 0.25
    ##   ..$ xGiven      : NULL
    ##   ..$ nExtra      : int 0
    ##   ..$ peakFinder  : NULL
    ##  $ sauve    :List of 7
    ##   ..$ minEval  : int 0
    ##   ..$ stateFile: NULL
    ##   ..$ flags    :List of 6
    ##   .. ..$ verbose   : int 0
    ##   .. ..$ final     : int 1
    ##   .. ..$ smooth    : int 0
    ##   .. ..$ keep_state: int 0
    ##   .. ..$ load_state: int 0
    ##   .. ..$ level     : int 0
    ##   ..$ rngSeed  : int 0
    ##   ..$ nNew     : int 1000
    ##   ..$ nMin     : int 50
    ##   ..$ flatness : num 50
    ##  $ vegas    :List of 8
    ##   ..$ minEval  : int 0
    ##   ..$ stateFile: NULL
    ##   ..$ flags    :List of 6
    ##   .. ..$ verbose   : int 0
    ##   .. ..$ final     : int 1
    ##   .. ..$ smooth    : int 0
    ##   .. ..$ keep_state: int 0
    ##   .. ..$ load_state: int 0
    ##   .. ..$ level     : int 0
    ##   ..$ rngSeed  : int 0
    ##   ..$ nStart   : int 1000
    ##   ..$ nIncrease: int 500
    ##   ..$ nBatch   : int 1000
    ##   ..$ gridNo   : int 0

## Vectorization

`cubintegrate` provides vector intefaces too: the parameter `nVec` is by
default 1, indicating a scalar interface. Any value \> 1 results in a
vectorized call. So `f` has to be constructed appropriately, thus:

``` r
my_dmvnorm_v <- function (x, mean, sigma, logdet) {
    distval <- stats::mahalanobis(t(x), center = mean, cov = sigma)
    exp(matrix(-(3 * log(2 * pi) + logdet + distval)/2, ncol = ncol(x)))
}
```

Here, the two underlying C libraries differ. The cubature library
manages the number of points used in vectorization dynamically and this
number can even vary from call to call. So any value of `nVec` greater
than 1 is merely a flag to use vectorization. The Cuba C library on the
other hand, will use the actual value of `nVec`.

``` r
cubintegrate(f = my_dmvnorm_v, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "pcubature",
             mean = rep(0, m), sigma = sigma, logdet = logdet,
             nVec = 128)
```

    ## $integral
    ## [1] 0.3341125
    ## 
    ## $error
    ## [1] 2.21207e-06
    ## 
    ## $neval
    ## [1] 4913
    ## 
    ## $returnCode
    ## [1] 0

``` r
cubintegrate(f = my_dmvnorm_v, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "cuhre",
             mean = rep(0, m), sigma = sigma, logdet = logdet,
             nVec = 128)
```

    ## $integral
    ## [1] 0.3341125
    ## 
    ## $error
    ## [1] 2.226266e-06
    ## 
    ## $nregions
    ## [1] 19
    ## 
    ## $neval
    ## [1] 4699
    ## 
    ## $prob
    ## [1] 0
    ## 
    ## $returnCode
    ## [1] 0

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
    ## [1] cubature_2.1.4-4
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    ##  [5] xfun_0.54         cachem_1.1.0      knitr_1.50        htmltools_0.5.8.1
    ##  [9] rmarkdown_2.30    lifecycle_1.0.4   cli_3.6.5         sass_0.4.10      
    ## [13] pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4   systemfonts_1.3.1
    ## [17] compiler_4.5.2    tools_4.5.2       ragg_1.5.0        bslib_0.9.0      
    ## [21] evaluate_1.0.5    Rcpp_1.1.0        yaml_2.3.10       jsonlite_2.0.0   
    ## [25] rlang_1.1.6       fs_1.6.6
