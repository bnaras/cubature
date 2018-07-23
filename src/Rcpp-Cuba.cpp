// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// Copyright (C) 2016 Balasubramanian Narasimhan
//
// We need both R and C/C++ interfaces!

#include <rcubature.h>
#include <cuba.h>

// Cuhre

// REMINDER nDim is 3, nComp = 1
int cuhre_fWrapper(const int *nDim, const double x[],
                  const int *nComp, double f[], void *userdata, const int *nVec,
                  const int *core) {

    //Rprintf("In Wrapper: nVec = %i\n", nVec);
    
    Rcpp::NumericVector xVal = Rcpp::NumericVector(x, x + (*nDim) * (*nVec));  /* The x argument for the R function f */
    //    Rcpp::Rcout<<"after xVal" <<std::endl;    
    if (*nVec > 1) {
        // Make the argument vector appear as a matrix for R
        xVal.attr("dim") = Rcpp::Dimension(*nDim, *nVec);
    }

    //Rcpp::Rcout<<"before call" <<std::endl;

    ii_ptr iip = (ii_ptr) userdata;
    Rcpp::NumericVector fx = Rcpp::Function(iip -> fun)(xVal);

    //Rcpp::Rcout<<"after call" <<std::endl;
    
    double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
    for (int i = 0; i < (*nComp) * (*nVec); ++i) {
        f[i] = fxp[i];
        // Rcpp::Rcout<< fval[i] <<std::endl;
    }
    // (iip -> count)++;
    return 0;
}

// int cuhre_fWrapper_v(const int *nDim, const double x[],
//                      const int *nComp, double f[], void *userdata, const int *nVec,
//                      const int *core) {

//     //Rcpp::Rcout << "In Wrapper: nVec = " << (*nVec) << std::endl;

//     Rcpp::NumericMatrix xVal(*nDim, *nVec, x);   /* The x argument for the R function f */

//     //Rcpp::Rcout<<"before call" <<std::endl;

//     ii_ptr iip = (ii_ptr) userdata;
//     Rcpp::NumericMatrix fx = Rcpp::Function(iip -> fun)(xVal);

//     //Rcpp::Rcout<<"after call" <<std::endl;

//     double* fxp = fx.begin();         /* The ptr to f(x) (real) matrix */
//     for (int i = 0; i < (*nComp) * (*nVec); ++i) {
//         f[i] = fxp[i];
//     }

//     return 0;
// }

// [[Rcpp::export]]
Rcpp::List doCuhre(int nComp, SEXP f, int nDim,
		   int nVec, int minEval, int maxEval, double absTol, double relTol,
                   int key, int flag) {

    Rcpp::NumericVector integral(nComp);
    Rcpp::NumericVector errVals(nComp);
    Rcpp::NumericVector prob(nComp);

    // Create a structure to hold integrand function and initialize it
    integrand_info II;
    II.fun = f;                 /* R function */

    int nregions, fail;

    // Set cores to be zero.
    cubacores(0, 0);
    
    //Rcpp::Rcout<<"Call Integrator" <<std::endl;
    //Rcpp::Rcout<<"Call Integrator nVec = " << nVec << std::endl;        
    Cuhre(nDim, nComp, (integrand_t) cuhre_fWrapper, (void *) &II, nVec,
          relTol, absTol, flag,
          minEval, maxEval, key,
          NULL, NULL,
          &nregions, &(II.count), &fail,
          integral.begin(), errVals.begin(), prob.begin());
    //Rcpp::Rcout<<"After Call Integrator" <<std::endl;

  return Rcpp::List::create(
			    Rcpp::_["integral"] = integral,
			    Rcpp::_["error"] = errVals,
			    Rcpp::_["nregions"] = nregions,
			    Rcpp::_["neval"] = II.count,
                            Rcpp::_["prob"] = prob,
			    Rcpp::_["returnCode"] = fail);
}

// Vegas

int vegas_or_suave_fWrapper(const int *nDim, const double x[],
                            const int *nComp, double f[], void *userdata, const int *nVec,
                            const int *core, const double weight[], const int *iter) {

    //    Rprintf("In Wrapper: nVec = %i\n", nVec);

    Rcpp::NumericVector xVal = Rcpp::NumericVector(x, x + (*nDim) * (*nVec));  /* The x argument for the R function f */
    //    Rcpp::Rcout<<"after xVal" <<std::endl;    
    if (*nVec > 1) {
        // Make the argument vector appear as a matrix for R
        xVal.attr("dim") = Rcpp::Dimension(*nDim, *nVec);
    }

    ii_ptr iip = (ii_ptr) userdata;
    Rcpp::NumericVector fx;

    //Rcpp::Rcout<<"before call" <<std::endl;        
    if (iip -> cuba_args) {
        Rcpp::NumericVector weightVal = Rcpp::NumericVector(weight, weight + (*nVec));  /* The weight argument for the R function f */
        //    Rcpp::Rcout<<"after weightVal" <<std::endl;    
        Rcpp::IntegerVector iterVal = Rcpp::IntegerVector(iter, iter + 1);  /* The iter argument for the R function f */        
        //    Rcpp::Rcout<<"before call" <<std::endl;
        fx = Rcpp::Function(iip -> fun)(xVal, Rcpp::_["cuba_weight"] = weightVal, Rcpp::_["cuba_iter"] = iterVal);
    } else {
        fx = Rcpp::Function(iip -> fun)(xVal);
    }

    //    Rcpp::Rcout<<"after call" <<std::endl;
    
    double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
    for (int i = 0; i < (*nComp) * (*nVec); ++i) {
        f[i] = fxp[i];
        //        Rcpp::Rcout<< f[i] <<std::endl;
    }
    return 0;
}

// int vegas_fWrapper_v(const int *nDim, const double x[],
//                      const int *nComp, double f[], void *userdata, const int *nVec,
//                      const int *core, const double weight[], const int *iter) {

//     //Rcpp::Rcout << "In Wrapper: nVec = " << (*nVec) << std::endl;

//     Rcpp::NumericMatrix xVal(*nDim, *nVec, x);   /* The x argument for the R function f */

//     ii_ptr iip = (ii_ptr) userdata;
//     Rcpp::NumericVector fx;

//     //Rcpp::Rcout<<"before call" <<std::endl;    
//     if (iip -> cuba_args) {
//         Rcpp::NumericVector weightVal = Rcpp::NumericVector(weight, weight + (*nVec));  /* The weight argument for the R function f */
//         Rcpp::IntegerVector iterVal = Rcpp::IntegerVector(iter, iter + 1);  /* The iter argument for the R function f */        
//         fx = Rcpp::Function(iip -> fun)(xVal, Rcpp::_["cuba_weight"] = weightVal, Rcpp::_["cuba_iter"] = iterVal);
//     } else {
//         fx = Rcpp::Function(iip -> fun)(xVal);
//     }
//     //Rcpp::Rcout<<"after call" <<std::endl;

//     double* fxp = fx.begin();         /* The ptr to f(x) (real) matrix */
//     for (int i = 0; i < (*nComp) * (*nVec); ++i) {
//         f[i] = fxp[i];
//     }

//     return 0;
// }

// [[Rcpp::export]]
Rcpp::List doVegas(int nComp, SEXP f, int nDim,
		   int nVec, int minEval, int maxEval, double absTol, double relTol,
                   int nStart, int nIncrease, int nBatch, int gridNo,
                   SEXP stateFile,
                   int seed, int flag, int cuba_args) {

    
    // Rcpp::Rcout<<"Entering" <<std::endl;
    Rcpp::NumericVector integral(nComp);
    Rcpp::NumericVector errVals(nComp);
    Rcpp::NumericVector prob(nComp);

    // Create a structure to hold integrand function and initialize it
    integrand_info II;
    II.cuba_args = cuba_args;   /* vegas specific args */
    II.fun = f;                 /* R function */

    int fail;

    // Set cores to be zero.
    cubacores(0, 0);
    
    char *filename = NULL;
    Rcpp::StringVector sv;
    if (!Rf_isNull(stateFile)) {
        sv = Rcpp::StringVector(stateFile);
        filename = sv(0);
        // Rcpp::Rcout << filename << std::endl;
    }
    //    Rcpp::Rcout<<"Call Integrator" <<std::endl;
    //        Rcpp::Rcout<<"Call Integrator nVec = " << nVec << std::endl;        
    Vegas(nDim, nComp, (integrand_t) vegas_or_suave_fWrapper, (void *) &II, nVec,
          relTol, absTol, flag,
          seed, minEval, maxEval,
          nStart, nIncrease, nBatch, gridNo, 
          filename,
          NULL,
          &(II.count), &fail,
          integral.begin(), errVals.begin(), prob.begin());
    //    Rcpp::Rcout<<"After Call Integrator" <<std::endl;

  return Rcpp::List::create(
			    Rcpp::_["integral"] = integral,
			    Rcpp::_["error"] = errVals,
			    Rcpp::_["neval"] = II.count,
                            Rcpp::_["prob"] = prob,
			    Rcpp::_["returnCode"] = fail);
}

// Suave

// int suave_fWrapper(const int *nDim, const double x[],
//                    const int *nComp, double f[], void *userdata, const int *nVec,
//                    const int *core, const double weight[], const int *iter) {

//     //    Rprintf("In Wrapper: nVec = %i\n", nVec);

//     Rcpp::NumericVector xVal = Rcpp::NumericVector(x, x + (*nDim) * (*nVec));  /* The x argument for the R function f */
//     //    Rcpp::Rcout<<"after xVal" <<std::endl;    
//     if (*nVec > 1) {
//         // Make the argument vector appear as a matrix for R
//         xVal.attr("dim") = Rcpp::Dimension(*nDim, *nVec);
//     }

//     ii_ptr iip = (ii_ptr) userdata;
//     Rcpp::NumericVector fx;

//     //    Rcpp::Rcout<<"before call" <<std::endl;        
//     if (iip -> cuba_args) {
//         Rcpp::NumericVector weightVal = Rcpp::NumericVector(weight, weight + (*nVec));  /* The weight argument for the R function f */
//         //    Rcpp::Rcout<<"after weightVal" <<std::endl;    
//         Rcpp::IntegerVector iterVal = Rcpp::IntegerVector(iter, iter + 1);  /* The iter argument for the R function f */        
//         //    Rcpp::Rcout<<"before call" <<std::endl;
//         fx = Rcpp::Function(iip -> fun)(xVal, Rcpp::_["cuba_weight"] = weightVal, Rcpp::_["cuba_iter"] = iterVal);
//     } else {
//         //        Rcpp::Rcout<<"No args" <<std::endl;                
//         fx = Rcpp::Function(iip -> fun)(xVal);
//     }

//     //    Rcpp::Rcout<<"after call" <<std::endl;
    
//     double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
//     for (int i = 0; i < (*nComp) * (*nVec); ++i) {
//         f[i] = fxp[i];
//         //        Rcpp::Rcout<< f[i] <<std::endl;
//     }
//     return 0;
// }

// int suave_fWrapper_v(const int *nDim, const double x[],
//                      const int *nComp, double f[], void *userdata, const int *nVec,
//                      const int *core, const double weight[], const int *iter) {

//     //Rcpp::Rcout << "In Wrapper: nVec = " << (*nVec) << std::endl;

//     Rcpp::NumericMatrix xVal(*nDim, *nVec, x);   /* The x argument for the R function f */

//     ii_ptr iip = (ii_ptr) userdata;
//     Rcpp::NumericVector fx;

//     //Rcpp::Rcout<<"before call" <<std::endl;    
//     if (iip -> cuba_args) {
//         Rcpp::NumericVector weightVal = Rcpp::NumericVector(weight, weight + (*nVec));  /* The weight argument for the R function f */
//         Rcpp::IntegerVector iterVal = Rcpp::IntegerVector(iter, iter + 1);  /* The iter argument for the R function f */        
//         fx = Rcpp::Function(iip -> fun)(xVal, Rcpp::_["cuba_weight"] = weightVal, Rcpp::_["cuba_iter"] = iterVal);
//     } else {
//         fx = Rcpp::Function(iip -> fun)(xVal);
//     }
//     //Rcpp::Rcout<<"after call" <<std::endl;

//     double* fxp = fx.begin();         /* The ptr to f(x) (real) matrix */
//     for (int i = 0; i < (*nComp) * (*nVec); ++i) {
//         f[i] = fxp[i];
//     }

//     return 0;
// }

// [[Rcpp::export]]
Rcpp::List doSuave(int nComp, SEXP f, int nDim,
		   int nVec, int minEval, int maxEval, double absTol, double relTol,
                   int nNew, int nMin, double flatness, 
                   SEXP stateFile,
                   int seed, int flag, int cuba_args) {

    
    //Rcpp::Rcout<<"Entering" <<std::endl;
    Rcpp::NumericVector integral(nComp);
    Rcpp::NumericVector errVals(nComp);
    Rcpp::NumericVector prob(nComp);

    // Create a structure to hold integrand function and initialize it
    integrand_info II;
    II.cuba_args = cuba_args; /* suave specific args */
    II.fun = f;                 /* R function */

    int nregions, fail;

    // Set cores to be zero.
    cubacores(0, 0);
    
    char *filename = NULL;
    Rcpp::StringVector sv;
    //Rcpp::Rcout << filename << std::endl;    
    if (!Rf_isNull(stateFile)) {
        sv = Rcpp::StringVector(stateFile);
        filename = sv(0);
        // Rcpp::Rcout << filename << std::endl;
    }
    //    Rcpp::Rcout<<"Call Integrator" <<std::endl;
    //        Rcpp::Rcout<<"Call Integrator nVec = " << nVec << std::endl;        
    Suave(nDim, nComp, (integrand_t) vegas_or_suave_fWrapper, (void *) &II, nVec,
          relTol, absTol, flag,
          seed, minEval, maxEval,
          nNew, nMin, flatness, 
          filename,
          NULL,
          &nregions, &(II.count), &fail,              
          integral.begin(), errVals.begin(), prob.begin());
    //    Rcpp::Rcout<<"After Call Integrator" <<std::endl;

  return Rcpp::List::create(
			    Rcpp::_["integral"] = integral,
			    Rcpp::_["error"] = errVals,
			    Rcpp::_["neval"] = II.count,
                            Rcpp::_["prob"] = prob,
			    Rcpp::_["returnCode"] = fail);
}

// Divonne

void peak_finder(const int *nDim, const double b[],
                 int *n, double x[], void *userdata) {

    ii_ptr iip = (ii_ptr) userdata;
        
    Rcpp::NumericMatrix bVal = Rcpp::NumericMatrix(2, *nDim, b); /* The b for the peakFinder */
    Rcpp::IntegerVector nVal = Rcpp::IntegerVector(n, n + 1);    /* The max number of pts to write */
    Rcpp::NumericMatrix fx = Rcpp::Function(iip -> peakFinder)(bVal, nVal);

    /* Update the actual number of points */
    *n = fx.nrow();
    
    double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
    for (int i = 0; i < (*n) * (*nDim); ++i) {
        x[i] = fxp[i];
    }
}

int divonne_fWrapper(const int *nDim, const double x[],
                   const int *nComp, double f[], void *userdata, const int *nVec,
                   const int *core, const int *phase) {

    //    Rprintf("In Wrapper: nVec = %i\n", nVec);

    Rcpp::NumericVector xVal = Rcpp::NumericVector(x, x + (*nDim) * (*nVec));  /* The x argument for the R function f */
    //    Rcpp::Rcout<<"after xVal" <<std::endl;    
    if (*nVec > 1) {
        // Make the argument vector appear as a matrix for R
        xVal.attr("dim") = Rcpp::Dimension(*nDim, *nVec);
    }

    ii_ptr iip = (ii_ptr) userdata;
    Rcpp::NumericVector fx;

    //    Rcpp::Rcout<<"before call" <<std::endl;        
    if (iip -> cuba_args) {
        Rcpp::IntegerVector phaseVal = Rcpp::IntegerVector(phase, phase + 1);  /* The phase argument for the R function f */        
        fx = Rcpp::Function(iip -> fun)(xVal, Rcpp::_["cuba_phase"] = phaseVal);
    } else {
        fx = Rcpp::Function(iip -> fun)(xVal);
    }
    double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
    for (int i = 0; i < (*nComp) * (*nVec); ++i) {
        f[i] = fxp[i];
    }
    return 0;
}

// [[Rcpp::export]]
Rcpp::List doDivonne(int nComp, SEXP f, int nDim,
                     int nVec, int minEval, int maxEval, double absTol, double relTol,
                     int key1, int key2, int key3,
                     int maxPass, double border,
                     double maxChisq, double minDeviation,
                     int nGiven, int ldxGiven,
                     SEXP xGiven, int nExtra,
                     SEXP peakFinder,
                     SEXP stateFile,
                     int seed, int flag, int cuba_args) {

    Rcpp::NumericVector integral(nComp);
    Rcpp::NumericVector errVals(nComp);
    Rcpp::NumericVector prob(nComp);

    // Create a structure to hold integrand function and initialize it
    integrand_info II;
    II.cuba_args = cuba_args;   /* cuba specific args */
    II.fun = f;                 /* R function */
    int peak_finder_null = Rf_isNull(peakFinder);
    if (!peak_finder_null) {
        II.fun = peakFinder;
    }
    
    int nregions, fail;

    // Set cores to be zero.
    cubacores(0, 0);
    
    char *filename = NULL;
    Rcpp::StringVector sv;
    //Rcpp::Rcout << filename << std::endl;    
    if (!Rf_isNull(stateFile)) {
        sv = Rcpp::StringVector(stateFile);
        filename = sv(0);
        // Rcpp::Rcout << filename << std::endl;
    }

    double* xGivenPtr = NULL; /* Ptr to xGiven start */
    Rcpp::NumericMatrix xGivenVal;
    if (!Rf_isNull(xGiven)) {
        xGivenVal = Rcpp::NumericMatrix(xGiven);
        xGivenPtr = xGivenVal.begin();
    }

    //    Rcpp::Rcout<<"Call Integrator" <<std::endl;
    Divonne(nDim, nComp, (integrand_t) divonne_fWrapper, (void *) &II, nVec,
            relTol, absTol, flag,
            seed, minEval, maxEval,
            key1, key2, key3,
            maxPass, border, maxChisq, minDeviation,
            nGiven, ldxGiven, xGivenPtr,
            nExtra,
            peak_finder_null ? NULL : (peakfinder_t) peak_finder,
            filename,
            NULL,
            &nregions, &(II.count), &fail,              
            integral.begin(), errVals.begin(), prob.begin());
    //    Rcpp::Rcout<<"After Call Integrator" <<std::endl;
    
    return Rcpp::List::create(
                              Rcpp::_["integral"] = integral,
                              Rcpp::_["error"] = errVals,
                              Rcpp::_["neval"] = II.count,
                              Rcpp::_["prob"] = prob,
                              Rcpp::_["returnCode"] = fail);
}




