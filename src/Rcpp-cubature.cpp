// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// Copyright (C) 2016 Balasubramanian Narasimhan
//
// We need both R and C/C++ interfaces!

#include <rcubature.h>
#include <cubature.h>

int fWrapper(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    //     Rcpp::Rcout<<"In Wrapper" <<std::endl;

    Rcpp::NumericVector xVal(x, x + ndim);   /* The x argument for the R function f */

    // Rcpp::Rcout<<"before call" <<std::endl;

    ii_ptr iip = (ii_ptr) fdata;
    Rcpp::NumericVector fx = Rcpp::Function(iip -> fun)(xVal);

    // Rcpp::Rcout<<"after call" <<std::endl;

    double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
    for (unsigned i = 0; i < fdim; ++i) {
        fval[i] = fxp[i];
    }
    (iip -> count)++;
    return 0;
}

int fWrapper_v(unsigned ndim, size_t npts, const double *x, void *fdata,
               unsigned fdim, double *fval) {
    //     Rcpp::Rcout<<"In Wrapper" <<std::endl;

    Rcpp::NumericMatrix xVal(ndim, npts, x);   /* The x argument for the R function f */

    //    Rcpp::Rcout<<"before call" <<std::endl;

    ii_ptr iip = (ii_ptr) fdata;
    Rcpp::NumericMatrix fx = Rcpp::Function(iip -> fun)(xVal);

    //    Rcpp::Rcout<<"after call" <<std::endl;

    double* fxp = fx.begin();         /* The ptr to f(x) (real) matrix */
    for (unsigned i = 0; i < fdim * npts; ++i) {
        fval[i] = fxp[i];
    }
    (iip -> count)++;
    return 0;
}

// [[Rcpp::export]]
Rcpp::List doHCubature(int fDim, SEXP f, Rcpp::NumericVector xLL, Rcpp::NumericVector xUL,
                       int maxEval, double absErr, double tol, int vectorInterface, unsigned norm) {

    Rcpp::NumericVector integral(fDim);
    Rcpp::NumericVector errVals(fDim);
    int retCode;

    // Create a structure to hold integrand function and initialize it
    integrand_info II;
    II.count = 0;               /* Zero count */    
    II.fun = f;                 /* R function */

    // Rcpp::Rcout<<"Call Integrator" <<std::endl;
    if (vectorInterface) {
        retCode = hcubature_v(fDim, fWrapper_v, (void *) &II,
                              xLL.size(), xLL.begin(), xUL.begin(),
                              maxEval, absErr, tol, (error_norm) norm,
                              integral.begin(), errVals.begin());
    } else {
        retCode = hcubature(fDim, fWrapper, (void *) &II,
                            xLL.size(), xLL.begin(), xUL.begin(),
                            maxEval, absErr, tol, (error_norm) norm,
                            integral.begin(), errVals.begin());
    }
    return Rcpp::List::create(
                              Rcpp::_["integral"] = integral,
                              Rcpp::_["error"] = errVals,
                              Rcpp::_["functionEvaluations"] = II.count,
                              Rcpp::_["returnCode"] = retCode);

}

// [[Rcpp::export]]
Rcpp::List doPCubature(int fDim, SEXP f, Rcpp::NumericVector xLL, Rcpp::NumericVector xUL,
                       int maxEval, double absErr, double tol, int vectorInterface, unsigned norm) {

    Rcpp::NumericVector integral(fDim);
    Rcpp::NumericVector errVals(fDim);
    int retCode;

    // Create a structure to hold integrand function and initialize it
    integrand_info II;
    II.count = 0;               /* Zero count */
    II.fun = f;                 /* R function */

    // Rcpp::Rcout<<"Call Integrator" <<std::endl;
    if (vectorInterface) {
        retCode = pcubature_v(fDim, fWrapper_v, (void *) &II,
                              xLL.size(), xLL.begin(), xUL.begin(),
                              maxEval, absErr, tol, (error_norm) norm,
                              integral.begin(), errVals.begin());
    } else {
        retCode = pcubature(fDim, fWrapper, (void *) &II,
                            xLL.size(), xLL.begin(), xUL.begin(),
                            maxEval, absErr, tol, (error_norm) norm,
                            integral.begin(), errVals.begin());
    }
    return Rcpp::List::create(
                              Rcpp::_["integral"] = integral,
                              Rcpp::_["error"] = errVals,
                              Rcpp::_["functionEvaluations"] = II.count,
                              Rcpp::_["returnCode"] = retCode);
}

