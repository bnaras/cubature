#' Experimental interface to newer Cuba library
#'
#' @importFrom Rcpp evalCpp
#'
#' @param f The function (integrand) to be integrated
#' @param lowerLimit The lower limit of integration, a vector for hypercubes
#' @param upperLimit The upper limit of integration, a vector for hypercubes
#' @param ...  All other arguments passed to the function f
#' @param tol The maximum tolerance, default 1e-5.
#' @param ncomp The dimension of the integrand, default 1, bears no relation to
#' the dimension of the hypercube
#' @param maxEval The maximum number of function evaluations needed, default 10^6.
#' Note that the actual number of function evaluations
#' performed is only approximately guaranteed not to exceed this number.
#' @param absError The maximum absolute error tolerated
#' @param key the quadrature rule key, default 0
#' @param nvec the number of vectorization points, ignored but fixed at 1 for now
#' @param doChecking A flag to be a bit anal about checking inputs to C
#' routines. A FALSE value results in approximately 9 percent speed gain in our
#' experiments. Your mileage will of course vary. Default value is FALSE.
#' @return The returned value is a list of five items:
#' \describe{
#'   \item{integral}{thevalue of the integral}
#'    \item{error}{the estimated relative error}
#'     \item{functionEvaluations}{the number of times the function was evaluated}
#'     \item{nregions}{the number of regions used}
#'     \item{returnCode}{the actual integer return code of the C routine}
#' }
#'
#' @export cuhre
#'
#' @author Balasubramanian Narasimhan
#' @references See \url{http://www.feynarts.de/cuba/}
#' @keywords math
cuhre <- function(f, lowerLimit, upperLimit, ..., tol = 1e-5,
                  ncomp = 1, maxEval = 10^6, absError = 0,
                  key = 0L, nvec = 1L,
                  doChecking = FALSE) {

    nL <- length(lowerLimit); nU <- length(upperLimit)
    if (ncomp <= 0 || nL <= 0 || nU <= 0) {
        stop("Both f and x must have dimension >= 1")
    }

    if (nL != nU) {
        stop("lowerLimit and upperLimit must have same length")
    }

    if (tol <= 0) {
        stop("tol should be positive!")
    }
    r <- upperLimit - lowerLimit

    f <- match.fun(f)

    if (doChecking) {
        fnF <- function(x) {
            fx <- f(lowerLimit + r * x, ...)
            if(!is.numeric(fx) || length(fx) != ncomp) {
                cat("hcubature: Error in evaluation function f(x) for x = ", x, "\n")
                stop("hcubature: Result f(x) is not numeric or has wrong dimension")
            }
            prod(r) * fx
        }
    } else {
        fnF <- function(x) {
            prod(r) * f(lowerLimit + r * x, ...)
        }
    }

    .Call('_cubature_doCuhre', PACKAGE = 'cubature', as.integer(ncomp), fnF,
          as.double(lowerLimit), as.double(upperLimit),
          as.integer(nvec), as.integer(maxEval),
          as.double(absError), as.double(tol), as.integer(key))
}
