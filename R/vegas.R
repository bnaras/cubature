#' Integration by a Monte Carlo Algorithm
#'
#' Implement a Monte Carlo algorithm for multidimensional numerical
#' integration.  This algorithm uses importance sampling as a
#' variance-reduction technique. Vegas iteratively builds up a
#' piecewise constant weight function, represented on a rectangular
#' grid. Each iteration consists of a sampling step followed by a
#' refinement of the grid.
#'
#' See details in the documentation.
#'
#' @importFrom Rcpp evalCpp
#'
#' @inheritParams cuhre
#' @param f The function (integrand) to be integrated as in \code{\link{cuhre}}. Optionally, the
#'     function can take two additional arguments in addition to the
#'     variable being integrated: - \code{cuba_weight} which is the
#'     weight of the point being sampled, - \code{cuba_iter} the
#'     current iteration number. The function author may choose to use
#'     these in any appropriate way or ignore them altogether.
#' @param flags flags governing the integration. A list with
#'     components: - \code{verbose}: \code{verbose} encodes the
#'     verbosity level, from 0 (default) to 3.  Level 0 does not print
#'     any output, level 1 prints \dQuote{reasonable} information on
#'     the progress of the integration, levels 2 and 3 echo the input
#'     parameters.  - \code{final}: when \code{ 0}, all sets of
#'     samples collected on a subregion during the various iterations
#'     or phases contribute to the final result.  When \code{ 1}, only
#'     the last (largest) set of samples is used in the final result.
#'     - \code{smooth}. When \code{smooth = 0}, apply additional
#'     smoothing to the importance function, this moderately improves
#'     convergence for many integrands.  When \code{smooth = 1} , use
#'     the importance function without smoothing, this should be
#'     chosen if the integrand has sharp edges.  - \code{keep_state}:
#'     when nonzero, retain state file if argument \code{stateFile} is
#'     non-null.  - \code{load_state}: when zero, load state file if
#'     found; if nonzero, reset state regardless - \code{level}: when
#'     \code{0}, Mersenne Twister random numbers are used. When
#'     nonzero Ranlux random numbers are used.  - \code{rngSeed}: When
#'     zero, Sobol quasi-random numbers are used for
#'     sampling. Otherwise the seed is used for the generator
#'     indicated by \code{level}.
#' @param nStart the number of integrand evaluations per iteration to
#'     start with.
#' @param nIncrease the increase in the number of integrand
#'     evaluations per iteration. The j-th iteration evaluates the
#'     integrand at nStart+(j-1)*nincrease points.
#' @param nBatch Vegas samples points not all at once, but in batches
#'     of a predetermined size, to avoid excessive memory
#'     consumption. \code{nbatch} is the number of points sampled in
#'     each batch. Tuning this number should usually not be necessary
#'     as performance is affected significantly only as far as the
#'     batch of samples fits into the CPU cache.
#' @param gridNo an integer.  Vegas may accelerate convergence to keep
#'     the grid accumulated during one integration for the next one,
#'     if the integrands are reasonably similar to each other. Vegas
#'     maintains an internal table with space for ten grids for this
#'     purpose.  If \code{gridno} is a number between 1 and 10, the
#'     grid is not discarded at the end of the integration, but stored
#'     in the respective slot of the table for a future
#'     invocation. The grid is only re-used if the dimension of the
#'     subsequent integration is the same as the one it originates
#'     from. In repeated invocations it may become necessary to flush
#'     a slot in memory. In this case the negative of the grid number
#'     should be set. Vegas will then start with a new grid and also
#'     restore the grid number to its positive value, such that at the
#'     end of the integration the grid is again stored in the
#'     indicated slot.
#' @seealso \code{\link{cuhre}}, \code{\link{suave}}, \code{\link{divonne}}
#'
#' @references G. P. Lepage (1978) A new algorithm for adaptive
#' multidimensional integration. \emph{J. Comput. Phys.}, \bold{27}, 192-210.
#'
#' G. P. Lepage (1980) VEGAS - An adaptive multi-dimensional integration
#' program. Research Report CLNS-80/447. Cornell University, Ithaca, N.-Y.
#'
#' T. Hahn (2005) CUBA-a library for multidimensional numerical integration.
#' \emph{Computer Physics Communications}, \bold{168}, 78-95.
#' @keywords math
#' @examples
#'
#' integrand <- function(arg, weight) {
#'   x <- arg[1]
#'   y <- arg[2]
#'   z <- arg[3]
#'   ff <- sin(x)*cos(y)*exp(z);
#' return(ff)
#' } # end integrand
#' vegas(integrand, lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
#'              relTol=1e-3,  absTol=1e-12,
#'              flags=list(verbose=2, final=0))
#'
#' @export vegas
vegas <- function(f, nComp = 1L, lowerLimit, upperLimit, ...,
                  relTol = 1e-5, absTol = 0,
                  minEval = 0L, maxEval = 10^6,
                  flags = list(verbose = 1, final = 1, smooth = 1, keep_state = 0, load_state = 0, level = 0, rngSeed = 12345L),
                  nVec = 1L, nStart = 1000L, nIncrease = 500L,
                  nBatch = 1000L, gridNo = 0L, stateFile = NULL) {

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
    all_flags <- list(verbose = 1, final = 1, smooth = 1, keep_state = 0, load_state = 0, level = 0, rngSeed = 12345L)
    for (x in names(flags)) all_flags[[x]] <- flags[[x]]

    f <- match.fun(f)
    cuba_params_exist <- length(intersect(c("cuba_weight", "cuba_iter"), names(formals(f)))) == 2L

    if (all(is.finite(c(lowerLimit, upperLimit)))) {
        r <- upperLimit - lowerLimit
        prodR <- prod(r)
        fnF <- function(x) prodR * f(lowerLimit + r * x, ...)
    } else {
        lowerLimit <- atan(lowerLimit)
        upperLimit <- atan(upperLimit)
        r <- upperLimit - lowerLimit
        prodR <- prod(r)
        fnF <- if (nVec > 1L)
                   if (cuba_params_exist) {
                       function(x, cuba_weight, cuba_iter) {
                           y <- lowerLimit + r * x
                           prodR * f(tan(y), cuba_weight = cuba_weight, cuba_iter = cuba_iter, ...) / rep(apply(cos(y), 2, prod)^2, each = nComp)
                       }
                   } else {
                       function(x) {
                           y <- lowerLimit + r * x
                           prodR * f(tan(y), ...) / rep(apply(cos(y), 2, prod)^2, each = nComp)
                       }
                   }
               else
                   if (cuba_params_exist) {
                       function(x, cuba_weight, cuba_iter) {
                       y <- lowerLimit + r * x
                       prodR * f(tan(y), cuba_weight = cuba_weight, cuba_iter = cuba_iter, ...) / prod(cos(y))^2
                       }
                   } else {
                       function(x) {
                           y <- lowerLimit + r * x
                           prodR * f(tan(y), ...) / prod(cos(y))^2
                       }
                   }
    }

    flag_code <- all_flags$verbose + 2^2 * all_flags$final + 2^3 * all_flags$smooth +
        2^4 * all_flags$keep_state + 2^5 * all_flags$load_state + 2^8 * all_flags$level

    .Call('_cubature_doVegas', PACKAGE = 'cubature',
          nComp, fnF, nL,
          nVec, minEval, maxEval,
          absTol, relTol, nStart,
          nIncrease, nBatch, gridNo,
          stateFile, all_flags$rngSeed, flag_code,
          cuba_params_exist)
}


