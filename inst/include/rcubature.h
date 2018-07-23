#include <Rcpp.h>      // need to include the main Rcpp header file only

typedef struct integrand_info {
  SEXP fun;                   /* The function itself */
  int count = 0;              /* Count of function evaluations */
  int cuba_args;              /* zero or 1 depending on cuba args present or not */
  SEXP peakFinder;            /* Peakfinder function for Divonne */
} *ii_ptr;


