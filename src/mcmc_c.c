#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>


// void F77_NAME(runMCMC_f)(int*N, double *result);
void F77_NAME(runMCMC_f)(int *cte, double *xvector, double *xvector_shifted);

extern SEXP c_runMCMC_f(SEXP cte, SEXP xvector){
  SEXP xvector_shifted;
  PROTECT(xvector_shifted = allocVector(REALSXP, 100));
  F77_CALL(runMCMC_f)(INTEGER(cte), REAL(xvector), REAL(xvector_shifted));
  UNPROTECT(1);
  return(xvector_shifted);
}

// extern SEXP c_runMCMC_f(SEXP N){
//   SEXP result;
//   PROTECT(result = allocVector(REALSXP, 100));
//   F77_CALL(runMCMC_f)(INTEGER(N), REAL(result));
//   UNPROTECT(1);
//   return(result);
// }

static const R_CallMethodDef CallEntries[] = {
  {"c_runMCMC_f",   (DL_FUNC) &c_runMCMC_f,  2},
  {NULL,         NULL,               0}
};

void R_init_dnIGEmodel (DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
