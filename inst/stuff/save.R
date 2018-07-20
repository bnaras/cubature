
#' @title Cuhre: Deterministic cubature rules
#'
#' @description Cuhre is a deterministic algorithm which uses one of
#'     several cubature rules of polynomial degree in a globally
#'     adaptive subdivision scheme. The subdivision algorithm is
#'     similar to [suave()] and works as follows:
#'
#' While the total estimated error exceeds the requested bounds:
#' 1. choose the region with the largest estimated error,
#' 2. bisect this region along the axis with the largest fourth difference,
#' 3. apply the cubature rule to the two subregions,
#' 4. merge the subregions into the list of regions and update the totals.
#'
#'Details on the algorithm and on the cubature rules employed in Cuhre can be found in the original references [@dcuhre]. The present
implementation offers only superficial improvements, such as an
interface consistent with the other [Cuba]{.smallcaps} routines and a
slightly simpler invocation, e.g. one does not have to allocate a
workspace.

In moderate dimensions Cuhre is very competitive, particularly if the
integrand is well approximated by polynomials. As the dimension
increases, the number of points sampled by the cubature rules rises
considerably, however, and by the same token the usefulness declines.
For the lower dimensions, the actual number of points that are spent per
invocation of the basic integration rule are listed in the following
table.


For vegas, it needs to take an additional argument named
#'     weight which is the weight of the point being sampled, which
#'     can either be ignored by the function or used in some
#'     constructive way. For divonne, the additional argument should
#'     be named phase and indicates the integration phase:
#' - 0, sampling of the points in \eqn{xgiven}:
#' - 1, partitioning phase,
#' - 2, final integration phase,
#' - 3, refinement phase.
#'     This information might be useful if the integrand takes long to
#'     compute and a sufficiently accurate approximation of the
#'     integrand is available. The actual value of the integral is
#'     only of minor importance in the partitioning phase, which is
#'     instead much more dependent on the peak structure of the
#'     integrand to find an appropriate tessellation. An approximation
#'     which reproduces the peak structure while leaving out the fine
#'     details might hence be a perfectly viable and much faster
#'     substitute when phase < 2. In all other instances, phase can be ignored.
