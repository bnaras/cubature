#' Integration by a Deterministic Iterative Adaptive Algorithm
#'
#' Implement a deterministic algorithm for multidimensional numerical
#' integration. Its algorithm uses one of several cubature rules in a globally
#' adaptive subdivision scheme.  The subdivision algorithm is similar to
#' \code{\link{suave}}'s.
#'
#' See details in the documentation.
#'
#' @param f The function (integrand) to be integrated. For cuhre, it
#'     can be something as simple as a function of a single argument,
#'     say x.
#' @param nComp The number of components of the integrand, default 1,
#'     bears no relation to the dimension of the hypercube over which
#'     integration is performed
#' @param lowerLimit The lower limit of integration, a vector for
#'     hypercubes
#' @param upperLimit The upper limit of integration, a vector for
#'     hypercubes
#' @param ...  All other arguments passed to the function f
#' @param relTol The maximum tolerance, default 1e-5.
#' @param absTol the absolute tolerance, default 0.
#' @param minEval the minimum number of function evaluations required
#' @param maxEval The maximum number of function evaluations needed,
#'     default 10^6.  Note that the actual number of function
#'     evaluations performed is only approximately guaranteed not to
#'     exceed this number.
#' @param key the quadrature rule key: \code{key = 7, 9, 11, 13}
#'     selects the cubature rule of degree key. Note that the
#'     degree-11 rule is available only in 3 dimensions, the degree-13
#'     rule only in 2 dimensions.
#' For other values, including the default 0, the rule is the degree-13
#' rule in 2 dimensions, the degree-11 rule in 3 dimensions, and the
#' degree-9 rule otherwise.
#' @param flags flags governing the integration. A list with components:
#' - \code{verbose}: \code{verbose} encodes the verbosity level, from 0 (default) to 3.
#' Level 0 does not print any output, level 1 prints \dQuote{reasonable}
#' information on the progress of the integration, level 2 also echoes the
#' input parameters, and level 3 further prints the subregion results.
#' - \code{final}: when \code{ 0}, all sets of samples collected on a subregion
#' during the various iterations or phases contribute to the final result.
#' When \code{ 1}, only the last (largest) set of samples is used in the final
#' result.
#' @param nVec the number of vectorization points, default 1, but can
#'     be set to an integer > 1 for vectorization, for example, 1024
#'     and the function f above needs to handle the vector of points
#'     appropriately
#' @param doChecking A flag to be a bit anal about checking inputs to
#'     C routines. A FALSE value can result in some performance gains
#'
#' @return A list with components:
#' \item{nregions }{the actual number of subregions needed}
#' \item{neval }{the actual number of integrand evaluations needed}
#' \item{returnCode}{return code: \code{0} , the desired accuracy was
#' reached, \code{-1}, dimension out of range, \code{1}, the accuracy
#' goal was not met within the allowed maximum number of integrand
#' evaluations.}  item{integral}{vector of length \code{nComp}; the
#' integral of \code{integrand} over the hypercube.}
#' \item{error}{vector of length \code{nComp}; the presumed absolute
#' error of \code{integral}} \item{prob}{vector of length
#' \code{nComp}; the \eqn{$\chi^2$}{Chi2}-probability (not the
#' \eqn{$\chi^2$}{Chi2}-value itself!) that \code{error} is not a
#' reliable estimate of the true integration error.}
#'
#' @seealso \code{\link{vegas}}, \code{\link{suave}}, \code{\link{divonne}}
#'
#' @importFrom Rcpp evalCpp
#'
#' @references J. Berntsen, T. O. Espelid (1991) An adaptive algorithm for the
#' approximate calculation of multiple integrals. \emph{ACM Transactions on
#' Mathematical Software}, \bold{17}(4), 437-451.
#'
#' T. Hahn (2005) CUBA-a library for multidimensional numerical integration.
#' \emph{Computer Physics Communications}, \bold{168}, 78-95.
#'
#' @references See \url{http://www.feynarts.de/cuba/}
#' @keywords math
#' @examples
#'
#' integrand <- function(arg) {
#'   x <- arg[1]
#'   y <- arg[2]
#'   z <- arg[3]
#'   ff <- sin(x)*cos(y)*exp(z);
#' return(ff)
#' } # End integrand
#'
#' NDIM <- 3
#' NCOMP <- 1
#' cuhre(f = integrand,
#'       lowerLimit = rep(0, NDIM),
#'       upperLimit = rep(1, NDIM),
#'       relTol = 1e-3, absTol= 1e-12,
#'       flags = list(verbose = 2, final = 0))
#'
#' @export cuhre
cuhre <- function(f, nComp = 1L, lowerLimit, upperLimit, ...,
                  relTol = 1e-5, absTol = 0,
                  minEval = 0L, maxEval = 10^6,
                  flags = list(verbose = 0, final = 1),
                  key = 0L, nVec = 1L,
                  doChecking = FALSE) {

    nL <- length(lowerLimit); nU <- length(upperLimit)
    if (nComp <= 0L || nL <= 0L || nU <= 0L) {
        stop("Both f and x must have dimension >= 1")
    }

    if (nL != nU) {
        stop("lowerLimit and upperLimit must have same length")
    }

    if (relTol <= 0) {
        stop("tol should be positive!")
    }
    r <- upperLimit - lowerLimit
    prodR <- prod(r)
    f <- match.fun(f)

    if (doChecking) {
        fnF <- function(x) {
            fx <- f(lowerLimit + r * x, ...)
            if(!is.numeric(fx) || length(fx) != nComp) {
                cat("cuhre: Error in evaluation function f(x) for x = ", x, "\n")
                stop("cuhre: Result f(x) is not numeric or has wrong dimension")
            }
            prodR * fx
        }
    } else {
        fnF <- function(x) {
            prodR  * f(lowerLimit + r * x, ...)
        }
    }

    flag_code <- flags$verbose + 4 * flags$final

    .Call('_cubature_doCuhre', PACKAGE = 'cubature',
          as.integer(nComp), fnF,
          as.double(lowerLimit), as.double(upperLimit),
          as.integer(nVec), as.integer(minEval), as.integer(maxEval),
          as.double(absTol), as.double(relTol), as.integer(key), as.integer(flag_code))
}


