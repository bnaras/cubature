#' Integration by Stratified Sampling for Variance Reduction
#'
#' Divonne works by stratified sampling, where the partioning of the
#' integration region is aided by methods from numerical optimization.
#'
#' Divonne uses stratified sampling for variance reduction, that is, it
#' partitions the integration region such that all subregions have an
#' approximately equal value of a quantity called the spread (volume times
#' half-range).
#'
#' See details in the documentation.
#'
#' @importFrom Rcpp evalCpp
#'
#' @inheritParams vegas
#' @param f The function (integrand) to be integrated as in
#'     \code{\link{cuhre}}. Optionally, the function can take an
#'     additional arguments in addition to the variable being
#'     integrated: - \code{cuba_phase}. The last argument, phase,
#'     indicates the integration phase: 0. sampling of the points in
#'     \code{xgiven}, 1. partitioning phase, 2. final integration
#'     phase, 3. refinement phase. This information might be useful if
#'     the integrand takes long to compute and a sufficiently accurate
#'     approximation of the integrand is available. The actual value
#'     of the integral is only of minor importance in the partitioning
#'     phase, which is instead much more dependent on the peak
#'     structure of the integrand to find an appropriate
#'     tessellation. An approximation which reproduces the peak
#'     structure while leaving out the fine details might hence be a
#'     perfectly viable and much faster substitute when
#'     \code{cuba_phase} < 2. In all other instances, phase can be
#'     ignored and it is entirely admissible to define the integrand
#'     without it. which is the
#' @param key1 integer that determines sampling in the partitioning
#'     phase: \code{key1 = 7, 9, 11, 13} selects the cubature rule of
#'     degree \code{key1}.  Note that the degree-11 rule is available
#'     only in 3 dimensions, the degree-13 rule only in 2
#'     dimensions. For other values of \code{key1}, a quasi-random
#'     sample of \eqn{n=|key1|}{\code{n=|key1|}} points is used, where
#'     the sign of \code{key1} determines the type of sample,
#'     \code{key1 = 0}, use the default rule. \code{key1 > 0}, use a
#'     Korobov quasi-random sample, \code{key1 < 0}, use a Sobol
#'     quasi-random sample if \code{flags$seed} is zero, otherwise a
#'     \dQuote{standard} sample (Mersenne Twister) pseudo-random
#'     sample
#' @param key2 integer that determines sampling in the final
#'     integration phase: same as \code{key1}, but here
#'     \eqn{$n=|key2|$}{\code{n = |key2|}} determines the number of
#'     points, \eqn{n > 39}{\code{n > 39}}, sample \eqn{n} points,
#'     \eqn{n < 40}{\code{n < 40}}, sample \eqn{n}{\code{n}}
#'     \code{nneed} points, where \code{nneed} is the number of points
#'     needed to reach the prescribed accuracy, as estimated by
#'     Divonne from the results of the partitioning phase.
#' @param key3 integer that sets the strategy for the refinement
#'     phase: \code{key3 = 0}, do not treat the subregion any further.
#'     \code{key3 = 1}, split the subregion up once more.  Otherwise,
#'     the subregion is sampled a third time with \code{key3}
#'     specifying the sampling parameters exactly as \code{key2}
#'     above.
#' @param maxPass integer that controls the thoroughness of the
#'     partitioning phase: The partitioning phase terminates when the
#'     estimated total number of integrand evaluations (partitioning
#'     plus final integration) does not decrease for \code{maxPass}
#'     successive iterations. A decrease in points generally indicates
#'     that Divonne discovered new structures of the integrand and was
#'     able to find a more effective partitioning. \code{maxPass} can
#'     be understood as the number of \dQuote{safety} iterations that
#'     are performed before the partition is accepted as final and
#'     counting consequently restarts at zero whenever new structures
#'     are found.
#' @param border the relative width of the border of the integration
#'     region.  Points falling into the border region will not be
#'     sampled directly, but will be extrapolated from two samples
#'     from the interior. Use a non-zero \code{border} if the
#'     integrand subroutine cannot produce values directly on the
#'     integration boundary. The relative width of the border is
#'     identical in all the dimensions. For example, set
#'     \code{border=0.1} for a border of width equal to 10\% of the
#'     width of the integration region.
#' @param maxChisq the maximum \eqn{\chi^2}{Chi2} value a single
#'     subregion is allowed to have in the final integration
#'     phase. Regions which fail this \eqn{\chi^2}{Chi2} test and
#'     whose sample averages differ by more than \code{min.deviation}
#'     move on to the refinement phase.
#' @param minDeviation a bound, given as the fraction of the requested
#'     error of the entire integral, which determines whether it is
#'     worthwhile further examining a region that failed the
#'     \eqn{\chi^2}{Chi2} test.  Only if the two sampling averages
#'     obtained for the region differ by more than this bound is the
#'     region further treated.
#' @param xGiven a matrix ( \code{nDim}, \code{nGiven}).  A list of
#'     \code{nGiven} points where the integrand might have peaks.
#'     Divonne will consider these points when partitioning the
#'     integration region.  The idea here is to help the integrator
#'     find the extrema of the integrand in the presence of very
#'     narrow peaks. Even if only the approximate location of such
#'     peaks is known, this can considerably speed up convergence.
#' @param nExtra the maximum number of extra points the peak-finder
#'     subroutine will return. If \code{nextra} is zero,
#'     \code{peakfinder} is not called and an arbitrary object may be
#'     passed in its place, e.g. just 0.
#' @param peakFinder the peak-finder subroutine. This R function is
#'     called whenever a region is up for subdivision and is supposed
#'     to point out possible peaks lying in the region, thus acting as
#'     the dynamic counterpart of the static list of points supplied
#'     in \code{xgiven}. It is expected to be declared as
#'     \code{peakFinder <- function(bounds, nMax)} where \code{bounds}
#'     is a matrix of dimension (\code{2, nDim}) which contains the
#'     lower (row 1) and upper (row 2) bounds of the subregion.  The
#'     returned value should be a matrix (\code{nX, nDim}) where
#'     \code{nX} is the actual number of points (should be less or
#'     equal to \code{nMax}).
#' @return A list with components: \describe{\item{nregions}{the actual
#'     number of subregions needed} \item{neval}{the actual number
#'     of integrand evaluations needed} \item{returnCode}{if zero,
#'     the desired accuracy was reached, if -1,
#'     dimension out of range, if 1, the accuracy goal was not met
#'     within the allowed maximum number of integrand evaluations.}
#'     \item{integral}{vector of length \code{nComp}; the integral of
#'     \code{integrand} over the hypercube} \item{error}{vector of
#'     length \code{nComp}; the presumed absolute error of
#'     \code{integral}} \item{prob}{vector of length \code{nComp};
#'     the \eqn{\chi^2}{Chi2}-probability (not the
#'     \eqn{\chi^2}{Chi2}-value itself!) that \code{error} is not a
#'     reliable estimate of the true integration error.}}
#'
#' @seealso \code{\link{cuhre}}, \code{\link{suave}}, \code{\link{vegas}}
#' @references J. H. Friedman, M. H. Wright (1981) A nested partitioning
#' procedure for numerical multiple integration. \emph{ACM Trans. Math.
#' Software}, \bold{7}(1), 76-92.
#'
#' J. H. Friedman, M. H. Wright (1981) User's guide for DIVONNE. SLAC Report
#' CGTM-193-REV, CGTM-193, Stanford University.
#'
#' T. Hahn (2005) CUBA-a library for multidimensional numerical integration.
#' \emph{Computer Physics Communications}, \bold{168}, 78-95.
#' @keywords math
#' @examples
#' integrand <- function(arg, phase) {
#'   x <- arg[1]
#'   y <- arg[2]
#'   z <- arg[3]
#'   ff <- sin(x)*cos(y)*exp(z);
#' return(ff)
#' }
#' divonne(integrand, relTol=1e-3,  absTol=1e-12, lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
#'         flags=list(verbose = 2),  key1= 47)
#'
#' # Example with a peak-finder function
#' nDim <- 3L
#' peakf <- function(bounds, nMax) {
#' #  print(bounds) # matrix (ndim,2)
#'   x <- matrix(0, ncol = nMax, nrow = nDim)
#'    pas <- 1 / (nMax - 1)
#'    # 1ier point
#'    x[, 1] <- rep(0, nDim)
#'    # Les autres points
#'    for (i in 2L:nMax) {
#'       x[, i] <- x[, (i - 1)] + pas
#'     }
#'   x
#' } #end peakf
#'
#' divonne(integrand, relTol=1e-3,  absTol=1e-12,
#'         lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
#'         flags=list(verbose = 2),  peakFinder = peakf, nExtra = 4L)
#' @export divonne
divonne <- function(f, nComp = 1L, lowerLimit, upperLimit, ...,
                    relTol = 1e-5, absTol = 1e-12,
                    minEval = 0L, maxEval = 10^6,
                    flags = list(verbose = 0L,
                                 final = 1L,
                                 keep_state = 0L,
                                 level = 0L),
                    rngSeed = 0L,
                    nVec = 1L,
                    key1 = 47L, key2 = 1L, key3 = 1L,
                    maxPass = 5L, border = 0, maxChisq = 10,
                    minDeviation = 0.25,
                    xGiven = NULL,  nExtra = 0L,
                    peakFinder = NULL,
                    stateFile = NULL) {
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

    ##  xGiven check
    if (!is.null(xGiven)) {
        if (!is.matrix(xGiven))
            stop("xGiven should be a matrix")
        nGiven <- ncol(xGiven)
        ldxGiven <- nrow(xGiven)
        if (ldxGiven != nL)
            stop("Matrix xgiven should have nDim rows")
    }
    else {
        nGiven <- 0
        ldxGiven <- nL
    }
    all_flags <- cuba_all_flags
    for (x in names(flags)) all_flags[[x]] <- flags[[x]]

    f <- match.fun(f)
    cuba_params_exist <- "cuba_phase" %in% names(formals(f))

    if (all(is.finite(c(lowerLimit, upperLimit)))) {
        r <- upperLimit - lowerLimit
        prodR <- prod(r)
        fnF <- if (cuba_params_exist) {
                   function(x, cuba_phase) prodR * f(lowerLimit + r * x, cuba_phase = cuba_phase, ...)
               } else {
                   function(x) prodR * f(lowerLimit + r * x, ...)
               }
        if (!is.null(xGiven)) {
            xGiven <- apply(xGiven, 1, function(x) x / r)
        }
    } else {
        lowerLimit <- atan(lowerLimit)
        upperLimit <- atan(upperLimit)
        r <- upperLimit - lowerLimit
        if (!is.null(xGiven)) {
            xGiven <- tan(xGiven)
        }
        prodR <- prod(r)
        fnF <- if (nVec > 1L) {
                   if (cuba_params_exist) {
                       function(x, cuba_phase) {
                           y <- lowerLimit + r * x
                           prodR * f(tan(y), cuba_phase = cuba_phase, ...) / rep(apply(cos(y), 2, prod)^2, each = nComp)
                       }
                   } else {
                       function(x) {
                           y <- lowerLimit + r * x
                           prodR * f(tan(y), ...) / rep(apply(cos(y), 2, prod)^2, each = nComp)
                       }
                   }
               } else {
                   if (cuba_params_exist) {
                       function(x, cuba_phase) {
                           y <- lowerLimit + r * x
                           prodR * f(tan(y), cuba_phase = cuba_phase, ...) / prod(cos(y))^2
                       }
                   } else {
                       function(x) {
                           y <- lowerLimit + r * x
                           prodR * f(tan(y), ...) / prod(cos(y))^2
                       }
                   }
               }
    }

    flag_code <- all_flags$verbose + 2^2 * all_flags$final +
        2^4 * all_flags$keep_state + 2^8 * all_flags$level

    .Call('_cubature_doDivonne', PACKAGE = 'cubature',
          nComp, fnF, nL,
          nVec, minEval, maxEval,
          absTol, relTol,
          key1, key2, key3,
          maxPass, border,
          maxChisq, minDeviation,
          nGiven, ldxGiven,
          xGiven, nExtra,
          peakFinder,
          stateFile,
          rngSeed, flag_code, cuba_params_exist)
}


