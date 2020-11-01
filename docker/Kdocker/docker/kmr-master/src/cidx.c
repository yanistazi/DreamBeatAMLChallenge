#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP cidx(SEXP x_, SEXP y_) {
  /*
   * C-Index for two real vectors x, y (UNCENSORED)
   * This is an alternative function that should output same number as Hmisc::rcorr.cens(x,y,outx=F)[[1]]
   */
  
  double *x, *y, s = 0;
  int i, j, n, n0 = 0;
  
  x = REAL(coerceVector(x_, REALSXP));
  y = REAL(coerceVector(y_, REALSXP));
  n = LENGTH(x_);
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (i != j && y[i] < y[j]) {
        ++n0;
        if (x[i] < x[j]) {
          s = s + 1;
        } else if (x[i] == x[j]) {
          s = s + 0.5;
        }
      }
    }
  }
  s = s/n0;
  
  return ScalarReal(s);
}
