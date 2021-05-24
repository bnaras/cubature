#ifndef __exp_cubature_h__
#define __exp_cubature_h__

#include <R.h>
#include <Rinternals.h>
#include "exp_cubature_typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* adapative integration by partitioning the integration domain ("h-adaptive")
     and using the same fixed-degree quadrature in each subdomain, recursively,
     until convergence is achieved. */
  int hcubature(unsigned fdim, integrand f, void *fdata,
		unsigned dim, const double *xmin, const double *xmax, 
		size_t maxEval, double reqAbsError, double reqRelError, 
		error_norm norm,
		double *val, double *err) {
    typedef int (*Fun)(unsigned, integrand, void*, unsigned, const double*,
		       const double*, size_t, double, double, error_norm,
		       double*, double*);
    static Fun fun = NULL;
    if (fun == NULL) {
      Rf_eval(Rf_lang2(Rf_install("loadNamespace"),
                       PROTECT(Rf_ScalarString(Rf_mkChar("cubature")))),
              R_GlobalEnv);
      UNPROTECT(1);
      fun = (Fun) R_GetCCallable("cubature", "hcubature");
    }
    return fun(fdim,f,fdata,dim,xmin,xmax,maxEval,reqAbsError,reqRelError,norm,val,err); 
  }
  
  /* as hcubature, but vectorized integrand */
  int hcubature_v(unsigned fdim, integrand_v f, void *fdata,
                  unsigned dim, const double *xmin, const double *xmax, 
                  size_t maxEval, double reqAbsError, double reqRelError, 
                  error_norm norm, double *val, double *err) {
    typedef int (*Fun)(unsigned, integrand_v, void*, unsigned, const double*,
                 const double*, size_t, double, double, error_norm,
                 double*, double*);
    static Fun fun = NULL;
    if (fun == NULL) {
      Rf_eval(Rf_lang2(Rf_install("loadNamespace"),
                       PROTECT(Rf_ScalarString(Rf_mkChar("cubature")))),
              R_GlobalEnv);
      UNPROTECT(1);
      fun = (Fun) R_GetCCallable("cubature", "hcubature_v");
    }
    return fun(fdim,f,fdata,dim,xmin,xmax,maxEval,reqAbsError,reqRelError,norm,val,err); 
  }
  
  /* adaptive integration by increasing the degree of (tensor-product
     Clenshaw-Curtis) quadrature rules ("p-adaptive"), rather than
     subdividing the domain ("h-adaptive").  Possibly better for
     smooth integrands in low dimensions. */
  int pcubature(unsigned fdim, integrand f, void *fdata,
                unsigned dim, const double *xmin, const double *xmax, 
                size_t maxEval, double reqAbsError, double reqRelError, 
                error_norm norm, double *val, double *err) {
    typedef int (*Fun)(unsigned, integrand, void*, unsigned, const double*,
                 const double*, size_t, double, double, error_norm,
                 double*, double*);
    static Fun fun = NULL;
    if (fun == NULL) {
      Rf_eval(Rf_lang2(Rf_install("loadNamespace"),
                       PROTECT(Rf_ScalarString(Rf_mkChar("cubature")))),
              R_GlobalEnv);
      UNPROTECT(1);
      fun = (Fun) R_GetCCallable("cubature", "pcubature");
    }
    return fun(fdim,f,fdata,dim,xmin,xmax,maxEval,reqAbsError,reqRelError,norm,val,err); 
  }
  
  /* as pcubature, but vectorized integrand */
  int pcubature_v(unsigned fdim, integrand_v f, void *fdata,
                  unsigned dim, const double *xmin, const double *xmax, 
                  size_t maxEval, double reqAbsError, double reqRelError, 
                  error_norm norm, double *val, double *err) {
    typedef int (*Fun)(unsigned, integrand_v, void*, unsigned, const double*,
                 const double*, size_t, double, double, error_norm,
                 double*, double*);
    static Fun fun = NULL;
    if (fun == NULL) {
      Rf_eval(Rf_lang2(Rf_install("loadNamespace"),
                       PROTECT(Rf_ScalarString(Rf_mkChar("cubature")))),
              R_GlobalEnv);
      UNPROTECT(1);
      fun = (Fun) R_GetCCallable("cubature", "pcubature_v");
    }
    return fun(fdim,f,fdata,dim,xmin,xmax,maxEval,reqAbsError,reqRelError,norm,val,err); 
  }
  
#ifdef __cplusplus
}
#endif

#endif
