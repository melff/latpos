#include "manifestos.h"

RcppExport SEXP ll_p (SEXP P_, SEXP Y_, SEXP N_, SEXP W_) {
BEGIN_RCPP
  using namespace Rcpp;
  
  NumericMatrix P(P_);
  NumericMatrix N(N_);
  NumericMatrix Y(Y_);
  NumericVector W(W_);

  int I = P.nrow();
  int J = P.ncol();

  if(N.nrow() != I || N.ncol() != J) return R_NilValue;
  if(W.size() != J) return R_NilValue;

  NumericMatrix Res(I,J);

  for(int i = 0; i < I; i++){

    Res(i,_) = W*Y(i,_)*N(i,_)*log(P(i,_));
  }
  return Res;
END_RCPP
}