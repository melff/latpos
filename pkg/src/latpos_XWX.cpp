#include "manifestos.h"

RcppExport SEXP latpos_XWX (SEXP X_, SEXP P_, SEXP N_, SEXP W_){
BEGIN_RCPP

  using namespace Rcpp;
  NumericMatrix X(X_), P(P_), N(N_);
  NumericVector W(W_);
  

  int J = P.ncol();
  int I = P.nrow();

  int K = X.ncol();

  arma::mat XWX = arma::zeros<arma::mat>(K,K);
  
  for(int j = 0; j < J; j++){

    int fromj = j*I;
    int toj = fromj + I-1;

    NumericMatrix Xj = X(Range(fromj,toj),_);
    NumericVector Pj = P(_,j);
    double wn = N(1,j)*W[j];
    
    arma::colvec Pvec = Rcpp::as<arma::colvec>(Pj);
    arma::mat WW = wn*(arma::diagmat(Pvec) - Pvec * trans(Pvec));
    arma::mat Xjmat = Rcpp::as<arma::mat>(Xj);

    XWX += trans(Xjmat) * WW * Xjmat;
  }
  return wrap(XWX);
END_RCPP
}

