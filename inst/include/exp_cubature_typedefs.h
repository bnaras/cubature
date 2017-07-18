#ifndef __exp_cubature_typedefs_h__
#define __exp_cubature_typedefs_h__

/* separate h file for the types and the error_norm enum,
 to avoid duplicate symbol error. */

#ifdef __cplusplus
extern "C" {
#endif
  
  /* a vector integrand - evaluates the function at the given point x
   (an array of length ndim) and returns the result in fval (an array
  of length fdim).   The void* parameter is there in case you have
  to pass any additional data through to your function (it corresponds
  to the fdata parameter you pass to cubature).  Return 0 on
  success or nonzero to terminate the integration. */
  typedef int (*integrand) (unsigned ndim, const double *x, void *,
               unsigned fdim, double *fval);
  
  /* a vector integrand of a vector of npt points: x[i*ndim + j] is the
  j-th coordinate of the i-th point, and the k-th function evaluation
  for the i-th point is returned in fval[i*fdim + k].  Return 0 on success
  or nonzero to terminate the integration. */
  typedef int (*integrand_v) (unsigned ndim, size_t npt,
               const double *x, void *,
               unsigned fdim, double *fval);
  
  /* Different ways of measuring the absolute and relative error when
  we have multiple integrands, given a vector e of error estimates
  in the individual components of a vector v of integrands.  These
  are all equivalent when there is only a single integrand. */
  typedef enum {
  ERROR_INDIVIDUAL = 0, /* individual relerr criteria in each component */
  ERROR_PAIRED, /* paired L2 norms of errors in each component,
    mainly for integrating vectors of complex numbers */
  ERROR_L2, /* abserr is L_2 norm |e|, and relerr is |e|/|v| */
  ERROR_L1, /* abserr is L_1 norm |e|, and relerr is |e|/|v| */
  ERROR_LINF /* abserr is L_\infty norm |e|, and relerr is |e|/|v| */
               } error_norm;
  
#ifdef __cplusplus
}
#endif

#endif
