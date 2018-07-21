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
#' @param f The function (integrand) to be integrated. Optionally, the
#'     function can take two additional arguments in addition to the
#'     variable being integrated: - \code{suave_weight} which is the
#'     weight of the point being sampled, - \code{suave_iter} the
#'     current iteration number. The function may choose to use these
#'     in any appropriate way or ignore them altogether
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
#' @param nNew the number of new integrand evaluations in each
#'     subdivision.
#' @param nMin the minimum number of samples a former pass must
#'     contribute to a subregion to be considered in that region’s
#'     compound integral value. Increasing nmin may reduce jumps in
#'     the \eqn{\Chi^2} value.
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
#' @param nVec the number of vectorization points, default 1, but can
#'     be set to an integer > 1 for vectorization, for example, 1024
#'     and the function f above needs to handle the vector of points
#'     appropriately
#' @param stateFile the name of an external file. Suave can store its
#'     entire internal state (i.e. all the information to resume an
#'     interrupted integration) in an external file.
#'
#' The state file is updated after every iteration. If, on a subsequent
#' invocation, Suave finds a file of the specified name, it loads the internal
#' state and continues from the point it left off. Needless to say, using an
#' existing state file with a different integrand generally leads to wrong
#' results. Once the integration finishes successfully, i.e. the prescribed
#' accuracy is attained, the state file is removed. This feature is useful
#' mainly to define \sQuote{check-points} in long-running integrations from
#' which the calculation can be restarted.
#' @return A list with components:
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
#' @seealso \code{\link{cuhre}}, \code{\link{divonne}}, \code{\link{vegas}}
#' @references T. Hahn (2005) CUBA-a library for multidimensional numerical
#' integration. \emph{Computer Physics Communications}, \bold{168}, 78-95.
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
#' suave(3, 1, integrand, rel.tol=1e-3,  abs.tol=1e-12,
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

    r <- upperLimit - lowerLimit
    prodR <- prod(r)

    f <- match.fun(f)
    dotargs_exist <- length(list(...)) > 0
    suave_params_exist <- length(intersect(c("suave_weight", "suave_iter"), names(formals(f)))) == 2L

    if (suave_params_exist && dotargs_exist) {
        fnF <- function(x, suave_weight, suave_iter) prodR * f(lowerLimit + r * x, suave_weight = suave_weight, suave_iter = suave_iter, ...)
    } else if (suave_params_exist && !dotargs_exist) {
        fnF <- function(x, suave_weight, suave_iter) prodR * f(lowerLimit + r * x, suave_weight = suave_weight, suave_iter = suave_iter)
    } else if (!suave_params_exist && dotargs_exist) {
        fnF <- function(x) prodR * f(lowerLimit + r * x, ...)
    } else {
        fnF <- function(x) prodR * f(lowerLimit + r * x)
    }

    flag_code <- all_flags$verbose + 2^2 * all_flags$final + 2^3 * all_flags$smooth +
        2^4 * all_flags$keep_state + 2^5 * all_flags$load_state + 2^8 * all_flags$level

    .Call('_cubature_doSuave', PACKAGE = 'cubature',
          nComp, fnF,
          lowerLimit, upperLimit,
          nVec, minEval, maxEval,
          absTol, relTol,
          nNew, nMin, flatness,
          stateFile, all_flags$rngSeed, flag_code,
          suave_params_exist)
}


