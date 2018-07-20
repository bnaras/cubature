// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// Copyright (C) 2016 Balasubramanian Narasimhan
//
// We need both R and C/C++ interfaces!

#include <rcubature.h>
#include <cuba.h>


// REMINDER nDim is 3, nComp = 1
int cuhre_fWrapper(const int *nDim, const double x[],
                  const int *nComp, double f[], void *userdata, const int *nVec,
                  const int *core) {

    //Rprintf("In Wrapper: nVec = %i\n", nVec);

    Rcpp::NumericVector xVal = Rcpp::NumericVector(x, x + (*nDim));  /* The x argument for the R function f */

    //Rcpp::Rcout<<"before call" <<std::endl;

    ii_ptr iip = (ii_ptr) userdata;
    Rcpp::NumericVector fx = Rcpp::Function(iip -> fun)(xVal);

    //Rcpp::Rcout<<"after call" <<std::endl;
    
    double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
    for (int i = 0; i < (*nComp); ++i) {
        f[i] = fxp[i];
        // Rcpp::Rcout<< fval[i] <<std::endl;
    }
    // (iip -> count)++;
    return 0;
}

int cuhre_fWrapper_v(const int *nDim, const double x[],
                     const int *nComp, double f[], void *userdata, const int *nVec,
                     const int *core) {

    //Rcpp::Rcout << "In Wrapper: nVec = " << (*nVec) << std::endl;

    Rcpp::NumericMatrix xVal(*nDim, *nVec, x);   /* The x argument for the R function f */

    //Rcpp::Rcout<<"before call" <<std::endl;

    ii_ptr iip = (ii_ptr) userdata;
    Rcpp::NumericMatrix fx = Rcpp::Function(iip -> fun)(xVal);

    //Rcpp::Rcout<<"after call" <<std::endl;

    double* fxp = fx.begin();         /* The ptr to f(x) (real) matrix */
    for (int i = 0; i < (*nComp) * (*nVec); ++i) {
        f[i] = fxp[i];
    }

    return 0;
}

// [[Rcpp::export]]
Rcpp::List doCuhre(int nComp, SEXP f, Rcpp::NumericVector xLL, Rcpp::NumericVector xUL,
		   int nVec, int minEval, int maxEval, double absTol, double relTol,
                   int key, int flag) {

    int nDim = xLL.size();
    Rcpp::NumericVector integral(nComp);
    Rcpp::NumericVector errVals(nComp);
    Rcpp::NumericVector prob(nComp);

    // Create a structure to hold integrand function and initialize it
    integrand_info II;
    // II.count = 0;               /* Zero count */
    II.fun = f;                 /* R function */

    int nregions, neval, fail;

    // Set cores to be zero.
    cubacores(0, 0);
    
    //Rcpp::Rcout<<"Call Integrator" <<std::endl;
    if (nVec > 1) {
        //Rcpp::Rcout<<"Call Integrator nVec = " << nVec << std::endl;        
        Cuhre(nDim, nComp, (integrand_t) cuhre_fWrapper_v, (void *) &II, nVec,
              relTol, absTol, flag,
              minEval, maxEval, key,
              NULL, NULL,
              &nregions, &neval, &fail,
              integral.begin(), errVals.begin(), prob.begin());
    } else {
        Cuhre(nDim, nComp, (integrand_t) cuhre_fWrapper, (void *) &II, nVec,
              relTol, absTol, flag,
              minEval, maxEval, key,
              NULL, NULL,
              &nregions, &neval, &fail,
              integral.begin(), errVals.begin(), prob.begin());
    }
    
    //Rcpp::Rcout<<"After Call Integrator" <<std::endl;

  return Rcpp::List::create(
			    Rcpp::_["integral"] = integral,
			    Rcpp::_["error"] = errVals,
			    Rcpp::_["nregions"] = nregions,
			    Rcpp::_["neval"] = neval,
                            Rcpp::_["prob"] = prob,
			    Rcpp::_["returnCode"] = fail);
}


int vegas_fWrapper(const int *nDim, const double x[],
                   const int *nComp, double f[], void *userdata, const int *nVec,
                   const int *core, const double weight[], const int *iter) {

    //    Rprintf("In Wrapper: nVec = %i\n", nVec);

    Rcpp::NumericVector xVal = Rcpp::NumericVector(x, x + (*nDim));  /* The x argument for the R function f */
    //    Rcpp::Rcout<<"after xVal" <<std::endl;
    Rcpp::NumericVector weightVal = Rcpp::NumericVector(weight, weight + (*nVec));  /* The weight argument for the R function f */
    //    Rcpp::Rcout<<"after weightVal" <<std::endl;    
    Rcpp::IntegerVector iterVal = Rcpp::IntegerVector(iter, iter + 1);  /* The iter argument for the R function f */        
    //    Rcpp::Rcout<<"before call" <<std::endl;

    ii_ptr iip = (ii_ptr) userdata;
    Rcpp::NumericVector fx = Rcpp::Function(iip -> fun)(xVal, Rcpp::_["vegas_weight"] = weightVal, Rcpp::_["vegas_iter"] = iterVal);

    //    Rcpp::NumericVector fx = Rcpp::Function(iip -> fun)(xVal);
    //    Rcpp::Rcout<<"after call" <<std::endl;
    
    double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
    for (int i = 0; i < (*nComp); ++i) {
        f[i] = fxp[i];
        //        Rcpp::Rcout<< f[i] <<std::endl;
    }
    return 0;
}

int vegas_fWrapper_v(const int *nDim, const double x[],
                     const int *nComp, double f[], void *userdata, const int *nVec,
                     const int *core, const double weight[], const int *iter) {

    //Rcpp::Rcout << "In Wrapper: nVec = " << (*nVec) << std::endl;

    Rcpp::NumericMatrix xVal(*nDim, *nVec, x);   /* The x argument for the R function f */
    Rcpp::NumericVector weightVal = Rcpp::NumericVector(weight, weight + (*nVec));  /* The weight argument for the R function f */
    Rcpp::IntegerVector iterVal = Rcpp::IntegerVector(iter, iter + 1);  /* The iter argument for the R function f */        

    //Rcpp::Rcout<<"before call" <<std::endl;

    ii_ptr iip = (ii_ptr) userdata;
    Rcpp::NumericVector fx = Rcpp::Function(iip -> fun)(xVal, weightVal, iterVal);    

    //Rcpp::Rcout<<"after call" <<std::endl;

    double* fxp = fx.begin();         /* The ptr to f(x) (real) matrix */
    for (int i = 0; i < (*nComp) * (*nVec); ++i) {
        f[i] = fxp[i];
    }

    return 0;
}

// [[Rcpp::export]]
Rcpp::List doVegas(int nComp, SEXP f, Rcpp::NumericVector xLL, Rcpp::NumericVector xUL,
		   int nVec, int minEval, int maxEval, double absTol, double relTol,
                   int nStart, int nIncrease, int nBatch, int gridNo, SEXP stateFile,
                   int seed, int flag) {

    
    Rcpp::Rcout<<"Entering" <<std::endl;
    int nDim = xLL.size();
    Rcpp::NumericVector integral(nComp);
    Rcpp::NumericVector errVals(nComp);
    Rcpp::NumericVector prob(nComp);

    // Create a structure to hold integrand function and initialize it
    integrand_info II;
    II.fun = f;                 /* R function */

    int neval, fail;

    // Set cores to be zero.
    cubacores(0, 0);
    
    //    Rcpp::Rcout<<"Call Integrator" <<std::endl;
    if (nVec > 1) {
        //        Rcpp::Rcout<<"Call Integrator nVec = " << nVec << std::endl;        
        Vegas(nDim, nComp, (integrand_t) vegas_fWrapper_v, (void *) &II, nVec,
              relTol, absTol, flag,
              seed, minEval, maxEval,
              nStart, nIncrease, nBatch, gridNo, 
              NULL, NULL,
              &neval, &fail,
              integral.begin(), errVals.begin(), prob.begin());
    } else {
        // Rcpp::Rcout<<"Call Integrator nVec = " << nVec << std::endl;        
        Vegas(nDim, nComp, (integrand_t) vegas_fWrapper, (void *) &II, nVec,
              relTol, absTol, flag,
              seed, minEval, maxEval,
              nStart, nIncrease, nBatch, gridNo, 
              NULL, NULL,
              &neval, &fail,
              integral.begin(), errVals.begin(), prob.begin());
    }
    
    //    Rcpp::Rcout<<"After Call Integrator" <<std::endl;

  return Rcpp::List::create(
			    Rcpp::_["integral"] = integral,
			    Rcpp::_["error"] = errVals,
			    Rcpp::_["neval"] = neval,
                            Rcpp::_["prob"] = prob,
			    Rcpp::_["returnCode"] = fail);
}
