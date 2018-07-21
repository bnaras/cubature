#include <Rcpp.h>      // need to include the main Rcpp header file only

typedef struct integrand_info {
    SEXP fun;                   /* The function itself */
    int count;                  /* Count of function evaluations, not used for Cuba */
    int vegas_args;             /* zero or 1 depending on vegas args present or not */
    int suave_args;             /* zero or 1 depending on suave args present or not */
    SEXP peakfinder;            /* Peakfinder function for Divonne */
} *ii_ptr;


