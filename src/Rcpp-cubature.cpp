// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// Copyright (C) 2016 Balasubramanian Narasimhan
//
// We need both R and C/C++ interfaces!

#include <Rcpp.h>      // need to include the main Rcpp header file only

#include "cubature.h"

SEXP fun;                   /* The function itself */
int count;                  /* Count of function evaluations */

void fWrapper(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    //     Rcpp::Rcout<<"In Wrapper" <<std::endl;

    Rcpp::NumericVector xVal(ndim);   /* The x argument for the R function f */
    double* xp = xVal.begin();        /* The ptr to x (real) vector */
    for (int i = 0; i < ndim; ++i) {
        xp[i] = x[i];
    }

    // Rcpp::Rcout<<"before call" <<std::endl;

    Rcpp::NumericVector fx = Rcpp::Function(fun)(xVal);

    // Rcpp::Rcout<<"after call" <<std::endl;

    double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
    for (int i = 0; i < fdim; ++i) {
        fval[i] = fxp[i];
    }
    count++;
}

// [[Rcpp::export]]
Rcpp::List doCubature(int fDim, SEXP f, Rcpp::NumericVector xLL, Rcpp::NumericVector xUL,
		      int maxEval, double absErr, double tol) {

    count = 0; /* Zero count */
    fun = f;

    Rcpp::NumericVector integral(fDim);
    Rcpp::NumericVector errVals(fDim);

    // Rcpp::Rcout<<"Call Integrator" <<std::endl;
    int retCode = adapt_integrate(fDim, fWrapper, NULL,
                                  xLL.size(), xLL.begin(), xUL.begin(),
                                  maxEval, absErr, tol,
                                  integral.begin(), errVals.begin());
    return Rcpp::List::create(
                              Rcpp::_["integral"] = integral,
                              Rcpp::_["error"] = errVals,
                              Rcpp::_["functionEvaluations"] = count,
                              Rcpp::_["returnCode"] = retCode);
}

