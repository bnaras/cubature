#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <cubature.h>


/* .Call entry points */

extern SEXP _cubature_doHCubature(SEXP fDimSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP maxEvalSEXP, SEXP absErrSEXP, SEXP tolSEXP, SEXP vectorInterfaceSEXP, SEXP normSEXP);

extern SEXP _cubature_doPCubature(SEXP fDimSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP maxEvalSEXP, SEXP absErrSEXP, SEXP tolSEXP, SEXP vectorInterfaceSEXP, SEXP normSEXP);


static const R_CallMethodDef CallEntries[] = {
  {"_cubature_doHCubature", (DL_FUNC) &_cubature_doHCubature, 9},
  {"_cubature_doPCubature", (DL_FUNC) &_cubature_doPCubature, 9},
  {NULL, NULL, 0}
};

void R_init_cubature(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);

  R_RegisterCCallable("cubature", "adapt_integrate", (DL_FUNC) hcubature);
  R_RegisterCCallable("cubature", "adapt_integrate_v", (DL_FUNC) hcubature_v);
  R_RegisterCCallable("cubature", "hcubature", (DL_FUNC) hcubature);
  R_RegisterCCallable("cubature", "hcubature_v", (DL_FUNC) hcubature_v);
  R_RegisterCCallable("cubature", "pcubature", (DL_FUNC) pcubature);
  R_RegisterCCallable("cubature", "pcubature_v", (DL_FUNC) pcubature_v);

}
