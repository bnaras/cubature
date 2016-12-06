#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>

#include "cubature.h"

/**
 * Register the native adapt_integrate functions so that they are C callable
 * (courtesy of Simen Gaure)
 */

void R_init_cubature(DllInfo *info) {
  R_RegisterCCallable("cubature", "adapt_integrate", (DL_FUNC) hcubature);
  R_RegisterCCallable("cubature", "adapt_integrate_v", (DL_FUNC) hcubature_v);
  R_RegisterCCallable("cubature", "hcubature", (DL_FUNC) hcubature);
  R_RegisterCCallable("cubature", "hcubature_v", (DL_FUNC) hcubature_v);
  R_RegisterCCallable("cubature", "pcubature", (DL_FUNC) pcubature);
  R_RegisterCCallable("cubature", "pcubature_v", (DL_FUNC) pcubature_v);

}
