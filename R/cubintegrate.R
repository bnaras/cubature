## Issue #17 by Simen Gaure
## Uniform interface for all methods.
##
cubature_common_defaults <- list(
    norm = c("INDIVIDUAL",
             "PAIRED",
             "L2",
             "L1",
             "LINF")
)

cuba_common_defaults <- list(
    minEval = 0L,
    stateFile = NULL
)

cuba_all_flags <- list(
    verbose = 0L,
    final = 1L,
    smooth = 0L,
    keep_state = 0L,
    load_state = 0L,
    level = 0L
)

#' Default arguments for each integration method
#'
#' Since each method has a different set of parameters, this function
#' returns the default values of all parameters that can be modified
#' and passed to integration routines.
#' @return a named list of parameters for each method.
#' @examples
#' default_args()
#' @export default_args
default_args <- function() {
    list(
        hcubature = cubature_common_defaults,
        pcubature = cubature_common_defaults,
        cuhre = c(cuba_common_defaults,
                  flags = list(cuba_all_flags),
                  list(
                      key = 0L
                  )
                  ),
        divonne = c(cuba_common_defaults,
                    flags = list(cuba_all_flags),
                    list(
                        rngSeed = 0L,
                        key1 = 47L,
                        key2 = 1L,
                        key3 = 1L,
                        maxPass = 5L,
                        border = 0,
                        maxChisq = 10,
                        minDeviation = 0.25,
                        xGiven = NULL,
                        nExtra = 0L,
                        peakFinder = NULL)
                    ),
        sauve = c(cuba_common_defaults,
                  flags = list(cuba_all_flags),
                  list(
                      rngSeed = 0L,
                      nNew = 1000L,
                      nMin = 50L,
                      flatness = 50)
                  ),
        vegas = c(cuba_common_defaults,
                  flags = list(cuba_all_flags),
                  list(
                      rngSeed = 0L,
                      nStart = 1000L,
                      nIncrease = 500L,
                      nBatch = 1000L,
                      gridNo = 0L)
                  )
    )
}

#' Unified Cubature Integration Interface
#'
#' Integrate a function within specified limits using method
#' specified. Further arguments specific to method as well as other
#' arguments to f may be passed. For defaults used in each method, see
#' help on the method or \code{\link{default_args}}'s.
#'
#' @param f The function (integrand) to be integrated. Can be
#'     vectorized version, but the additional arguments \code{...}
#'     must indicate via either \code{vectorInterface = TRUE} for
#'     `hcubature` and `pcubature`, or a value for \code{nVec}. See
#'     details on each method.
#' @param fDim The number of components of f, default 1, bears no
#'     relation to the dimension of the hypercube over which
#'     integration is performed.
#' @param lower The lower limit of integration, a vector for
#'     hypercubes.
#' @param upper The upper limit of integration, a vector for
#'     hypercubes.
#' @param relTol The maximum tolerance, default 1e-5.
#' @param absTol the absolute tolerance, default 1e-12.
#' @param maxEval The maximum number of function evaluations needed,
#'     default 10^6.  Note that the actual number of function
#'     evaluations performed is only approximately guaranteed not to
#'     exceed this number.
#' @param nVec the number of vectorization points for Cuba C library,
#'     default 1, but can be set to an integer > 1 for vectorization,
#'     for example, 1024. The function f above needs to handle the
#'     vector of points appropriately; see vignette examples. Unlike
#'     Cuba, the cubature C library manages the number of points on
#'     its own and can vary between calls. Therefore, any value for
#'     nVec greater than one implies vectorization for a cubature
#'     method.
#' @param method the method to use should be one of "hcubature",
#'     "pcubature", "cuhre", "divonne", "suave" or "vegas".
#' @param ...  All other arguments which may include integration
#'     method specific parameters and those for f. Unrecognized
#'     parameters for integration method are presumed to be intended
#'     for f and so processed.
#' @return The returned value is a list of items: -\item{integral}{the
#'     value of the integral} - \item{error}{the estimated relative
#'     error for cubature; for Cuba it is the estimated absolute
#'     error} \item{neval}{the number of times the function was
#'     evaluated} - \item{returnCode}{the actual integer return code
#'     of the C routine; a non-zero value usually indicates problems;
#'     further interpretation depends on method} - \item{nregions}{for
#'     Cuba routines, the actual number of subregions needed}
#'     \item{prob}{the \eqn{\chi^2}{Chi2}-probability (not the
#'     \eqn{\chi^2}{Chi2}-value itself!) that \code{error} is not a
#'     reliable estimate of the true integration error.}
#'
#' @seealso \code{\link{default_args}}, \code{\link{hcubature}},
#'     \code{\link{pcubature}}, \code{\link{cuhre}},
#'     \code{\link{vegas}}, \code{\link{suave}}, \code{\link{divonne}}
#'
#' @examples
#' I.1d <- function(x) {
#'   sin(4*x) *
#'     x * ((x * ( x * (x*x-4) + 1) - 1))
#' }
#' I.1d_v <- function(x) {
#'    matrix(apply(x, 2, function(z)
#'        sin(4 * z) *
#'        z * ((z * ( z * (z * z - 4) + 1) - 1))),
#'        ncol = ncol(x))
#' }
#' cubintegrate(f = I.1d, lower = -2, upper = 2, method = "pcubature")
#' cubintegrate(f = I.1d, lower = -2, upper = 2, method = "cuhre", flags=list(verbose = 2))
#' cubintegrate(f = I.1d_v, lower = -2, upper = 2, method = "hcubature", nVec = 2L)
#' cubintegrate(f = I.1d_v, lower = -2, upper = 2, method = "cuhre", nVec = 128L)
#'
#' @export cubintegrate
cubintegrate <- function(f, lower, upper, fDim = 1,
                         method = c('hcubature',
                                    'pcubature',
                                    'cuhre',
                                    'divonne',
                                    'suave',
                                    'vegas'),
                         relTol = 1e-5,
                         absTol = 1e-12,
                         maxEval = 10^6,
                         nVec = 1L,
                         ...) {
    method <- match.arg(method)
    other_args <- list(...)
    result <-
        if (grepl("cubature", method)) {
            do.call(method,
                    args = c(list(f = f, fDim = fDim, lowerLimit = lower, upperLimit = upper),
                             list(tol = relTol, absError = absTol, maxEval = maxEval,
                                  vectorInterface = (nVec > 1L)),
                             other_args))

        } else {
            do.call(method,
                    args = c(list(f = f, nComp = fDim, lowerLimit = lower, upperLimit = upper),
                             list(relTol = relTol, absTol = absTol, maxEval = maxEval, nVec = nVec),
                             other_args))
        }

    names(result) <- gsub("functionEvaluations", "neval", names(result))
    result
}
