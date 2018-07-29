#' Integration with SUbregion-Adaptive Vegas Algorithm
#'
#' Suave uses \code{\link{vegas}}-like importance sampling combined with a
#' globally adaptive subdivision strategy: Until the requested accuracy is
#' reached, the region with the largest error at the time is bisected in the
#' dimension in which the fluctuations of the integrand are reduced most. The
#' number of new samples in each half is prorated for the fluctuation in that
#' half.
#'
#' See details in the documentation.
#'
#' @importFrom Rcpp evalCpp
#'
#' @inheritParams vegas
#' @param nNew the number of new integrand evaluations in each
#'     subdivision.
#' @param nMin the minimum number of samples a former pass must
#'     contribute to a subregion to be considered in that region’s
#'     compound integral value. Increasing nmin may reduce jumps in
#'     the \eqn{\chi^2}{Chi2} value.
#' @param flatness the parameter p, or the type of norm used to
#'     compute the fluctuation of a sample. This determines how
#'     prominently ‘outliers,’ i.e. individual samples with a large
#'     fluctuation, figure in the total fluctuation, which in turn
#'     determines how a region is split up. As suggested by its name,
#'     flatness should be chosen large for ‘flat’ integrands and small
#'     for ‘volatile’ integrands with high peaks. Note that since
#'     flatness appears in the exponent, one should not use too large
#'     values (say, no more than a few hundred) lest terms be
#'     truncated internally to prevent overflow.
#'
#' @seealso \code{\link{cuhre}}, \code{\link{divonne}}, \code{\link{vegas}}
#' @references T. Hahn (2005) CUBA-a library for multidimensional numerical
#' integration. \emph{Computer Physics Communications}, \bold{168}, 78-95.
#' @keywords math
#' @examples
#'
#' integrand <- function(arg) {
#'   x <- arg[1]
#'   y <- arg[2]
#'   z <- arg[3]
#'   ff <- sin(x)*cos(y)*exp(z);
#' return(ff)
#' } # end integrand
#' suave(integrand, lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
#'              relTol=1e-3,  absTol=1e-12,
#'              flags=list(verbose=2, final=0))
#'
#' @export suave
suave <- function(f, nComp = 1L, lowerLimit, upperLimit, ...,
                  relTol = 1e-5, absTol = 0,
                  minEval = 0L, maxEval = 10^6,
                  flags = list(verbose = 1, final = 1, smooth = 1, keep_state = 0, load_state = 0, level = 0, rngSeed = 12345L),
                  nVec = 1L, nNew = 1000L, nMin = 50L,
                  flatness = 50, stateFile = NULL) {

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

    flag_code <- all_flags$verbose + 2^2 * all_flags$final + 2^3 * all_flags$smooth +
        2^4 * all_flags$keep_state + 2^5 * all_flags$load_state + 2^8 * all_flags$level

    .Call('_cubature_doSuave', PACKAGE = 'cubature',
          nComp, fnF, nL,
          nVec, minEval, maxEval,
          absTol, relTol,
          nNew, nMin, flatness,
          stateFile, all_flags$rngSeed, flag_code,
          cuba_params_exist)
}


