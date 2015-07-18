#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>

#include "cubature.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("mylib", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif


SEXP CUB_common_env;   /* The common environment we use */
SEXP f;                /* The function itself */
int count;             /* Count of function evaluations */

SEXP CUB_set_common_env(SEXP rho) {
    if (!isEnvironment(rho))
        error(_("Argument rho must be an environment"));
    CUB_common_env = rho;
    return R_NilValue;
}

void fWrapper(unsigned ndim, const double *x, void *fdata,
	      unsigned fdim, double *fval) {
  SEXP xx, fx;
  double *rx, *rfx;
  int i;

  PROTECT(xx = allocVector(REALSXP, ndim));
  rx = REAL(xx);
  for (i = 0; i < ndim; ++i) {
    rx[i] = x[i];
  }
  defineVar(install("x"), xx, CUB_common_env);
  PROTECT(fx = eval(f, CUB_common_env));

  rfx = REAL(fx);
  for (i = 0; i < fdim; ++i) {
    fval[i] = rfx[i];
  }
  UNPROTECT(2);
  count++;
}

SEXP doCubature(SEXP sfDim, SEXP sf, SEXP sxLL, SEXP sxUL, SEXP smaxEval,
		SEXP sabsErr, SEXP stol, SEXP rho) {

  double *xLL, *xUL, *val, *err;
  double absErr, tol;
  int i, fDim, nDim, maxEval, retCode;
  SEXP integral, errVals, fCount, rCode, ans;

  /* Save the environment for later use */
  CUB_common_env = rho;
  f = sf;
  count = 0; /* zero the count of function evaluations */

  fDim = INTEGER(sfDim)[0]; nDim = LENGTH(sxLL);
  xLL = REAL(sxLL); xUL = REAL(sxUL);
  absErr = REAL(sabsErr)[0]; tol = REAL(stol)[0];
  maxEval = INTEGER(smaxEval)[0];

  val = (double *) R_alloc(fDim, sizeof(double));
  err = (double *) R_alloc(fDim, sizeof(double));

  retCode = adapt_integrate(fDim, fWrapper, NULL,
			    nDim, xLL, xUL, maxEval,
			    absErr, tol,
			    val, err);

  PROTECT(integral = allocVector(REALSXP, fDim));
  for (i = 0; i < fDim; ++i) {
    REAL(integral)[i] = val[i];
  }

  PROTECT(errVals = allocVector(REALSXP, fDim));
  for (i = 0; i < fDim; ++i) {
    REAL(errVals)[i] = err[i];
  }

  PROTECT(fCount = allocVector(INTSXP, 1));
  INTEGER(fCount)[0] = count;

  PROTECT(rCode = allocVector(INTSXP, 1));
  INTEGER(rCode)[0] = retCode;

  PROTECT(ans = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ans, 0, integral);
  SET_VECTOR_ELT(ans, 1, errVals);
  SET_VECTOR_ELT(ans, 2, fCount);
  SET_VECTOR_ELT(ans, 3, rCode);

  UNPROTECT(5);
  return ans;
}

/**
 * Register the native adapt_integrate functions so that they are C callable
 * (courtesy of Simen Gaure)
 */

void R_init_cubature(DllInfo *info) {
  R_RegisterCCallable("cubature", "adapt_integrate", (DL_FUNC) adapt_integrate);
  R_RegisterCCallable("cubature", "adapt_integrate_v", (DL_FUNC) adapt_integrate_v);
}
