#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <cubature.h>


/* .Call entry points */

extern SEXP _cubature_doHCubature(SEXP fDimSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP maxEvalSEXP, SEXP absErrSEXP, SEXP tolSEXP, SEXP vectorInterfaceSEXP, SEXP normSEXP);

extern SEXP _cubature_doPCubature(SEXP fDimSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP maxEvalSEXP, SEXP absErrSEXP, SEXP tolSEXP, SEXP vectorInterfaceSEXP, SEXP normSEXP);

extern SEXP _cubature_doCuhre(SEXP nCompSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP nVecSEXP, SEXP minEvalSEXP, SEXP maxEvalSEXP, SEXP absTolSEXP, SEXP relTolSEXP, SEXP keySEXP, SEXP flagSEXP);

extern SEXP _cubature_doVegas(SEXP nCompSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP nVecSEXP, SEXP minEvalSEXP, SEXP maxEvalSEXP, SEXP absTolSEXP, SEXP relTolSEXP, SEXP nStartSEXP, SEXP nIncreaseSEXP, SEXP nBatchSEXP, SEXP gridNoSEXP, SEXP stateFileSEXP, SEXP seedSEXP, SEXP flagSEXP, SEXP vegas_argsSEXP);

extern SEXP _cubature_doSuave(SEXP nCompSEXP, SEXP fSEXP, SEXP xLLSEXP, SEXP xULSEXP, SEXP nVecSEXP, SEXP minEvalSEXP, SEXP maxEvalSEXP, SEXP absTolSEXP, SEXP relTolSEXP, SEXP nNewSEXP, SEXP nMinSEXP, SEXP flatnessSEXP, SEXP stateFileSEXP, SEXP seedSEXP, SEXP flagSEXP, SEXP suave_argsSEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_cubature_doHCubature", (DL_FUNC) &_cubature_doHCubature, 9},
  {"_cubature_doPCubature", (DL_FUNC) &_cubature_doPCubature, 9},
  {"_cubature_doCuhre", (DL_FUNC) &_cubature_doCuhre, 11},
  {"_cubature_doVegas", (DL_FUNC) &_cubature_doVegas, 17},
  {"_cubature_doSuave", (DL_FUNC) &_cubature_doSuave, 16},  
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
