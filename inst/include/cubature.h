/* Adaptive multidimensional integration of a vector of integrands.
 *
 * Copyright (c) 2005-2013 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 *     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
 *     (http://www.gnu.org/software/gsl/)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef CUBATURE_H
#define CUBATURE_H

#include <stdlib.h> /* for size_t */

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* USAGE: Call hcubature or pcubature with your function as described
          in the README file. */

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

/* Integrate the function f from xmin[dim] to xmax[dim], with at most
   maxEval function evaluations (0 for no limit), until the given
   absolute or relative error is achieved.  val returns the integral,
   and err returns the estimate for the absolute error in val; both
   of these are arrays of length fdim, the dimension of the vector
   integrand f(x). The return value of the function is one of:

     CUBATURE_SUCCESS         (0) - converged to the requested tolerance
     CUBATURE_FAILURE         (1) - internal error (e.g. allocation failure,
                                    invalid parameters); val/err are unset
     CUBATURE_NOT_CONVERGED   (2) - the maxEval budget was exhausted before
                                    the requested tolerance was reached;
                                    val/err hold the best available estimate
                                    and its (possibly large) error bound,
                                    which the caller should treat as
                                    best-effort, not guaranteed

   Note: prior versions returned CUBATURE_SUCCESS even when maxEval was
   exhausted without meeting the tolerance. Callers that relied on that
   silent behavior should now explicitly check for CUBATURE_NOT_CONVERGED
   and act accordingly. */
#define CUBATURE_SUCCESS        0
#define CUBATURE_FAILURE        1
#define CUBATURE_NOT_CONVERGED  2

/* adapative integration by partitioning the integration domain ("h-adaptive")
   and using the same fixed-degree quadrature in each subdomain, recursively,
   until convergence is achieved. */
int hcubature(unsigned fdim, integrand f, void *fdata,
	      unsigned dim, const double *xmin, const double *xmax,
	      size_t maxEval, double reqAbsError, double reqRelError,
	      error_norm norm,
	      double *val, double *err);

/* as hcubature, but vectorized integrand */
int hcubature_v(unsigned fdim, integrand_v f, void *fdata,
		unsigned dim, const double *xmin, const double *xmax,
		size_t maxEval, double reqAbsError, double reqRelError,
		error_norm norm,
		double *val, double *err);

/* as hcubature / hcubature_v, but with optional "robust" error
   estimation enabled. When robust != 0, three additional safeguards
   are active, at the cost of some extra flops per region and a
   modest increase in function evaluations for non-smooth integrands:
   (1) a degree-3 null rule difference is computed alongside the
       existing degree-5 one and combined with a Cuhre-style
       decay-check formula (Berntsen, Espelid & Genz, ACM TOMS 17(4),
       437-451 (1991));
   (2) a parent-child consistency check is applied at the split step,
       inflating the children's error by the discrepancy between the
       sum of their estimates and the parent's prior estimate;
   (3) a suspect-region detector based on running peak density flags
       large regions whose density is far below the observed peak and
       forces further subdivision.
   When robust == 0 the behavior is identical to hcubature /
   hcubature_v, bit-for-bit, and no extra work is done.
   These safeguards address the TODO at the top of hcubature.c (two-
   level error estimation) and substantially reduce "fool's
   convergence" on integrands with thin or localized support. */
int hcubature_robust(unsigned fdim, integrand f, void *fdata,
		     unsigned dim, const double *xmin, const double *xmax,
		     size_t maxEval, double reqAbsError, double reqRelError,
		     error_norm norm,
		     double *val, double *err, int robust);
int hcubature_v_robust(unsigned fdim, integrand_v f, void *fdata,
		       unsigned dim, const double *xmin, const double *xmax,
		       size_t maxEval, double reqAbsError, double reqRelError,
		       error_norm norm,
		       double *val, double *err, int robust);

/* adaptive integration by increasing the degree of (tensor-product
   Clenshaw-Curtis) quadrature rules ("p-adaptive"), rather than
   subdividing the domain ("h-adaptive").  Possibly better for
   smooth integrands in low dimensions. */
int pcubature_v_buf(unsigned fdim, integrand_v f, void *fdata,
		    unsigned dim, const double *xmin, const double *xmax,
		    size_t maxEval, 
		    double reqAbsError, double reqRelError,
		    error_norm norm,
		    unsigned *m,
		    double **buf, size_t *nbuf, size_t max_nbuf,
		    double *val, double *err);
int pcubature_v(unsigned fdim, integrand_v f, void *fdata,
		unsigned dim, const double *xmin, const double *xmax, 
		size_t maxEval, double reqAbsError, double reqRelError, 
		error_norm norm,
		double *val, double *err);
int pcubature(unsigned fdim, integrand f, void *fdata,
	      unsigned dim, const double *xmin, const double *xmax,
	      size_t maxEval, double reqAbsError, double reqRelError,
	      error_norm norm,
	      double *val, double *err);

/* Robust variants of pcubature / pcubature_v: when robust > 0, the
   Clenshaw-Curtis rule starts at order m=robust in every dimension
   (instead of m=0), providing denser initial sampling that is less
   likely to structurally miss localized integrand features. Cost
   scales as (2^(robust+1)+1)^dim extra function evaluations up
   front. robust = 0 is identical to the non-robust entry points.
   This does NOT promise correctness for arbitrarily localized
   integrands; for those, prefer hcubature(..., robust = 1) or
   Cuba cuhre / divonne. */
int pcubature_robust(unsigned fdim, integrand f, void *fdata,
		     unsigned dim, const double *xmin, const double *xmax,
		     size_t maxEval, double reqAbsError, double reqRelError,
		     error_norm norm,
		     double *val, double *err, int robust);
int pcubature_v_robust(unsigned fdim, integrand_v f, void *fdata,
		       unsigned dim, const double *xmin, const double *xmax,
		       size_t maxEval, double reqAbsError, double reqRelError,
		       error_norm norm,
		       double *val, double *err, int robust);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* CUBATURE_H */
