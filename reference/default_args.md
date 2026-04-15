# Default arguments for each integration method

Since each method has a different set of parameters, this function
returns the default values of all parameters that can be modified and
passed to integration routines.

## Usage

``` r
default_args()
```

## Value

a named list of parameters for each method.

## Examples

``` r
default_args()
#> $hcubature
#> $hcubature$norm
#> [1] "INDIVIDUAL" "PAIRED"     "L2"         "L1"         "LINF"      
#> 
#> 
#> $pcubature
#> $pcubature$norm
#> [1] "INDIVIDUAL" "PAIRED"     "L2"         "L1"         "LINF"      
#> 
#> 
#> $cuhre
#> $cuhre$minEval
#> [1] 0
#> 
#> $cuhre$stateFile
#> NULL
#> 
#> $cuhre$flags
#> $cuhre$flags$verbose
#> [1] 0
#> 
#> $cuhre$flags$final
#> [1] 1
#> 
#> $cuhre$flags$smooth
#> [1] 0
#> 
#> $cuhre$flags$keep_state
#> [1] 0
#> 
#> $cuhre$flags$load_state
#> [1] 0
#> 
#> $cuhre$flags$level
#> [1] 0
#> 
#> 
#> $cuhre$key
#> [1] 0
#> 
#> 
#> $divonne
#> $divonne$minEval
#> [1] 0
#> 
#> $divonne$stateFile
#> NULL
#> 
#> $divonne$flags
#> $divonne$flags$verbose
#> [1] 0
#> 
#> $divonne$flags$final
#> [1] 1
#> 
#> $divonne$flags$smooth
#> [1] 0
#> 
#> $divonne$flags$keep_state
#> [1] 0
#> 
#> $divonne$flags$load_state
#> [1] 0
#> 
#> $divonne$flags$level
#> [1] 0
#> 
#> 
#> $divonne$rngSeed
#> [1] 0
#> 
#> $divonne$key1
#> [1] 47
#> 
#> $divonne$key2
#> [1] 1
#> 
#> $divonne$key3
#> [1] 1
#> 
#> $divonne$maxPass
#> [1] 5
#> 
#> $divonne$border
#> [1] 0
#> 
#> $divonne$maxChisq
#> [1] 10
#> 
#> $divonne$minDeviation
#> [1] 0.25
#> 
#> $divonne$xGiven
#> NULL
#> 
#> $divonne$nExtra
#> [1] 0
#> 
#> $divonne$peakFinder
#> NULL
#> 
#> 
#> $suave
#> $suave$minEval
#> [1] 0
#> 
#> $suave$stateFile
#> NULL
#> 
#> $suave$flags
#> $suave$flags$verbose
#> [1] 0
#> 
#> $suave$flags$final
#> [1] 1
#> 
#> $suave$flags$smooth
#> [1] 0
#> 
#> $suave$flags$keep_state
#> [1] 0
#> 
#> $suave$flags$load_state
#> [1] 0
#> 
#> $suave$flags$level
#> [1] 0
#> 
#> 
#> $suave$rngSeed
#> [1] 0
#> 
#> $suave$nNew
#> [1] 1000
#> 
#> $suave$nMin
#> [1] 50
#> 
#> $suave$flatness
#> [1] 50
#> 
#> 
#> $vegas
#> $vegas$minEval
#> [1] 0
#> 
#> $vegas$stateFile
#> NULL
#> 
#> $vegas$flags
#> $vegas$flags$verbose
#> [1] 0
#> 
#> $vegas$flags$final
#> [1] 1
#> 
#> $vegas$flags$smooth
#> [1] 0
#> 
#> $vegas$flags$keep_state
#> [1] 0
#> 
#> $vegas$flags$load_state
#> [1] 0
#> 
#> $vegas$flags$level
#> [1] 0
#> 
#> 
#> $vegas$rngSeed
#> [1] 0
#> 
#> $vegas$nStart
#> [1] 1000
#> 
#> $vegas$nIncrease
#> [1] 500
#> 
#> $vegas$nBatch
#> [1] 1000
#> 
#> $vegas$gridNo
#> [1] 0
#> 
#> 
```
