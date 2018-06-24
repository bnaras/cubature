// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// Copyright (C) 2016 Balasubramanian Narasimhan
//
// We need both R and C/C++ interfaces!

#include <rcubature.h>
#include <cuba.h>


// REMINDER ndim is 3, ncomp = 1
int cuba_fWrapper(const int *ndim, const double x[],
                  const int *ncomp, double f[], void *userdata, const int *nvec,
                  const int *spin) {

    //Rprintf("In Wrapper: nvec = %i\n", nvec);

    Rcpp::NumericVector xVal(*ndim);   /* The x argument for the R function f */
    double* xp = xVal.begin();        /* The ptr to x (real) vector */
    for (int i = 0; i < (*ndim); ++i) {
        xp[i] = x[i];
        // Rcpp::Rcout<< x[j][i] <<std::endl;
    }

    //Rcpp::Rcout<<"before call" <<std::endl;

    ii_ptr iip = (ii_ptr) userdata;
    Rcpp::NumericVector fx = Rcpp::Function(iip -> fun)(xVal);

    //Rcpp::Rcout<<"after call" <<std::endl;
    
    double* fxp = fx.begin();         /* The ptr to f(x) (real) vector */
    for (int i = 0; i < (*ncomp); ++i) {
        f[i] = fxp[i];
        // Rcpp::Rcout<< fval[i] <<std::endl;
    }
    (iip -> count)++;
    return 0;
}

int cuba_fWrapper_v(const int *ndim, const double x[],
                    const int *ncomp, double f[], void *userdata, const int *nvec) {

    //Rcpp::Rcout << "In Wrapper: nvec = " << (*nvec) << std::endl;

    Rcpp::NumericMatrix xVal(*ndim, *nvec);   /* The x argument for the R function f */
    double* xp = xVal.begin();        /* The ptr to x (real) matrix */
    for (int i = 0; i < (*ndim) * (*nvec); ++i) {
        xp[i] = x[i];
    }

    //Rcpp::Rcout<<"before call" <<std::endl;

    ii_ptr iip = (ii_ptr) userdata;
    Rcpp::NumericMatrix fx = Rcpp::Function(iip -> fun)(xVal);

    //Rcpp::Rcout<<"after call" <<std::endl;

    double* fxp = fx.begin();         /* The ptr to f(x) (real) matrix */
    for (int i = 0; i < (*ncomp) * (*nvec); ++i) {
        f[i] = fxp[i];
    }

    return 0;
}

// [[Rcpp::export]]
Rcpp::List doCuhre(int ncomp, SEXP f, Rcpp::NumericVector xLL, Rcpp::NumericVector xUL,
		   int nvec, int maxEval, double absErr, double tol, int key) {

    int ndim = xLL.size();
    Rcpp::NumericVector integral(ncomp);
    Rcpp::NumericVector errVals(ncomp);
    Rcpp::NumericVector prob(ncomp);

    // Create a structure to hold integrand function and initialize it
    integrand_info II;
    II.count = 0;               /* Zero count */
    II.fun = f;                 /* R function */

    int nregions, neval, fail;

    // Set cores to be zero.
    cubacores(0, 0);
    
    //Rcpp::Rcout<<"Call Integrator" <<std::endl;
    if (nvec > 1) {
        //Rcpp::Rcout<<"Call Integrator nvec = " << nvec << std::endl;        
        Cuhre(ndim, ncomp, (integrand_t) cuba_fWrapper_v, (void *) &II, nvec,
              tol, absErr, 0,
              0, maxEval, key,
              NULL, NULL,
              &nregions, &neval, &fail,
              integral.begin(), errVals.begin(), prob.begin());
    } else {
        Cuhre(ndim, ncomp, (integrand_t) cuba_fWrapper, (void *) &II, nvec,
              tol, absErr, 0,
              0, maxEval, key,
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
			    Rcpp::_["returnCode"] = fail);
}

