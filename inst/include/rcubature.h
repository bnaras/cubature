#include <Rcpp.h>      // need to include the main Rcpp header file only

typedef struct integrand_info {
    SEXP fun;                   /* The function itself */
    int count;                  /* Count of function evaluations */
} *ii_ptr;


