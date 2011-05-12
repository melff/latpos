#include "manifestos.h"

RcppExport SEXP latpos_p (SEXP A_, SEXP B_){

  using namespace Rcpp;
  NumericMatrix A(A_), B(B_);

  if (A.ncol() != B.ncol()) return R_NilValue;

  int m = A.nrow();
  int n = B.nrow();
  int ndims = A.ncol();

  NumericMatrix P(m,n);

  for(int i = 0; i < m; i++){

    for(int j = 0; j < n; j++){

      double tmp = 0.0;
      for(int d = 0; d < ndims; d++){

        double diff = A(i,d) - B(j,d);
        tmp += diff*diff/2;
      }
      P(i,j) = exp(-tmp);
    }
  }

  for(int j = 0; j < n; j++){

    double tmp = 0.0;
    for(int i = 0; i < m; i++){


      tmp += P(i,j);
    }
    for(int i = 0; i < m; i++){

      P(i,j) = P(i,j)/tmp;
    }

  }

  return P;
}
