#include "manifestos.h"

RcppExport SEXP latpos_resid (SEXP P_, SEXP Y_, SEXP N_, SEXP W_) {
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

  NumericVector Res(I*J);
  
  for(int j = 0; j < J; j++){

    Range r(j*I,j*I+I-1);
    
    Res[r] = (W[j]*N(_,j))*(Y(_,j)-P(_,j));
  }
  return Res;
END_RCPP
}