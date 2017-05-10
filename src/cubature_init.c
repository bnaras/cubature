#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


/* .Call entry points */

extern SEXP cubature_doHCubature(SEXP fDimSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP maxEvalSEXP, SEXP absErrSEXP, SEXP tolSEXP, SEXP vectorInterfaceSEXP, SEXP normSEXP);

extern SEXP cubature_doPCubature(SEXP fDimSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP maxEvalSEXP, SEXP absErrSEXP, SEXP tolSEXP, SEXP vectorInterfaceSEXP, SEXP normSEXP);


static const R_CallMethodDef CallEntries[] = {
  {"cubature_doHCubature", (DL_FUNC) &cubature_doHCubature, 9},
  {"cubature_doPCubature", (DL_FUNC) &cubature_doPCubature, 9},
  {NULL, NULL, 0}
};

void R_init_cubature(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
