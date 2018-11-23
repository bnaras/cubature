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
#' @param nComp The number of components of f, default 1, bears no
#'     relation to the dimension of the hypercube over which
#'     integration is performed.
#' @param lowerLimit The lower limit of integration, a vector for
#'     hypercubes.
#' @param upperLimit The upper limit of integration, a vector for
#'     hypercubes.
#' @param ...  All other arguments passed to the function f.
#' @param relTol The maximum tolerance, default 1e-5.
#' @param absTol the absolute tolerance, default 1e-12.
#' @param minEval the minimum number of function evaluations required
#' @param maxEval The maximum number of function evaluations needed,
#'     default 10^6.  Note that the actual number of function
#'     evaluations performed is only approximately guaranteed not to
#'     exceed this number.
#' @param key the quadrature rule key: \code{key = 7, 9, 11, 13}
#'     selects the cubature rule of degree key. Note that the
#'     degree-11 rule is available only in 3 dimensions, the degree-13
#'     rule only in 2 dimensions.  For other values, including the
#'     default 0, the rule is the degree-13 rule in 2 dimensions, the
#'     degree-11 rule in 3 dimensions, and the degree-9 rule
#'     otherwise.
#' @param flags flags governing the integration. The list here is
#'     exhaustive to keep the documentation and invocation uniform,
#'     but not all flags may be used for a particular method as noted
#'     below.  List components: \describe{ \item{\code{verbose}}{
#'     encodes the verbosity level, from 0 (default) to 3.  Level 0
#'     does not print any output, level 1 prints reasonable
#'     information on the progress of the integration, level 2 also
#'     echoes the input parameters, and level 3 further prints the
#'     subregion results.}  \item{\code{final}}{when 0, all sets of
#'     samples collected on a subregion during the various iterations
#'     or phases contribute to the final result.  When 1, only the
#'     last (largest) set of samples is used in the final result.}
#'     \item{\code{smooth}}{Applies to Suave and Vegas only. When 0,
#'     apply additional smoothing to the importance function, this
#'     moderately improves convergence for many integrands.  When 1,
#'     use the importance function without smoothing, this should be
#'     chosen if the integrand has sharp edges.}
#'     \item{\code{keep_state}}{when nonzero, retain state file if
#'     argument \code{stateFile} is non-null, else delete
#'     \code{stateFile} if specified.}
#'     \item{\code{load_state}}{Applies to Vegas only. Reset the
#'     integrator's state even if a state file is present, i.e. keep
#'     only the grid. Together with \code{keep_state} this allows a
#'     grid adapted by one integration to be used for another
#'     integrand.}  \item{\code{level}}{applies only to Divonne, Suave
#'     and Vegas. When \code{0}, Mersenne Twister random numbers are
#'     used. When nonzero Ranlux random numbers are used, except when
#'     \code{rngSeed} is zero which forces use of Sobol quasi-random
#'     numbers. Ranlux implements Marsaglia and Zaman's 24-bit RCARRY
#'     algorithm with generation period p, i.e. for every 24 generated
#'     numbers used, another p-24 are skipped. The luxury level for
#'     the Ranlux generator may be encoded in \code{level} as follows:
#'     \describe{ \item{Level 1 (p = 48)}{gives very long period,
#'     passes the gap test but fails spectral test} \item{Level 2 (p =
#'     97)}{passes all known tests, but theoretically still defective}
#'     \item{Level 3 (p = 223)}{any theoretically possible
#'     correlations have very small chance of being observed}
#'     \item{Level 4 (p = 389)}{highest possible luxury, all 24 bits
#'     chaotic} \item{Levels 5-23}{default to 3, values above 24
#'     directly specify the period p.}}  Note that Ranlux's original
#'     level 0, (mis)used for selecting Mersenne Twister in Cuba, is
#'     equivalent to \code{level} = 24.}}
#' @param nVec the number of vectorization points, default 1, but can
#'     be set to an integer > 1 for vectorization, for example, 1024
#'     and the function f above needs to handle the vector of points
#'     appropriately. See vignette examples.
#' @param stateFile the name of an external file. Vegas can store its
#'     entire internal state (i.e. all the information to resume an
#'     interrupted integration) in an external file.  The state file
#'     is updated after every iteration. If, on a subsequent
#'     invocation, Vegas finds a file of the specified name, it loads
#'     the internal state and continues from the point it left
#'     off. Needless to say, using an existing state file with a
#'     different integrand generally leads to wrong results. Once the
#'     integration finishes successfully, i.e. the prescribed accuracy
#'     is attained, the state file is removed. This feature is useful
#'     mainly to define \sQuote{check-points} in long-running
#'     integrations from which the calculation can be restarted.
#' @return A list with components: \describe{\item{nregions}{the
#'     actual number of subregions needed} \item{neval}{the actual
#'     number of integrand evaluations needed} \item{returnCode}{if
#'     zero, the desired accuracy was reached, if -1, dimension out of
#'     range, if 1, the accuracy goal was not met within the allowed
#'     maximum number of integrand evaluations.}
#'     \item{integral}{vector of length \code{nComp}; the integral of
#'     \code{integrand} over the hypercube} \item{error}{vector of
#'     length \code{nComp}; the presumed absolute error of
#'     \code{integral}} \item{prob}{vector of length \code{nComp}; the
#'     \eqn{\chi^2}{Chi2}-probability (not the
#'     \eqn{\chi^2}{Chi2}-value itself!) that \code{error} is not a
#'     reliable estimate of the true integration error.}}
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
                  relTol = 1e-5, absTol = 1e-12,
                  minEval = 0L, maxEval = 10^6,
                  flags = list(verbose = 0L,
                               final = 1L,
                               keep_state = 0L,
                               level = 0L),
                  key = 0L, nVec = 1L, stateFile = NULL) {

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
    f <- match.fun(f)

    all_flags <- cuba_all_flags
    for (x in names(flags)) all_flags[[x]] <- flags[[x]]

    if (all(is.finite(c(lowerLimit, upperLimit)))) {
        r <- upperLimit - lowerLimit
        prodR <- prod(r)
        fnF <- function(x) {
            prodR  * f(lowerLimit + r * x, ...)
        }
    } else {
        lowerLimit <- atan(lowerLimit)
        upperLimit <- atan(upperLimit)
        r <- upperLimit - lowerLimit
        prodR <- prod(r)
        fnF <- if (nVec > 1L)
                   function(x) {
                       y <- lowerLimit + r * x
                       prodR * f(tan(y), ...) / rep(apply(cos(y), 2, prod)^2, each = nComp)
                   }
        else
                   function(x) {
                       y <- lowerLimit + r * x
                       prodR * f(tan(y), ...) / prod(cos(y))^2
                   }
    }

    flag_code <- all_flags$verbose + 2^2 * all_flags$final +
        2^4 * all_flags$keep_state

    .Call('_cubature_doCuhre', PACKAGE = 'cubature',
          nComp, fnF, nL,
          nVec, minEval, maxEval,
          absTol, relTol,
          stateFile,
          key, flag_code)
}


