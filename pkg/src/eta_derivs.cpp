#include "manifestos.h"

RcppExport SEXP d_eta_d_phi (SEXP A_, SEXP B_, SEXP Q_){

  using namespace Rcpp;
  NumericMatrix A(A_), B(B_), Q(Q_);

  if (A.ncol() != B.ncol()) return R_NilValue;

  int m = A.nrow();
  int n = B.nrow();
  int ndims = A.ncol();
  int m_ndims = m*ndims;

  if(Q.nrow()!=m_ndims) return R_NilValue;
  int r = Q.ncol();
  
  NumericMatrix X(n*m,r);
  NumericMatrix Tmp(m,ndims);
  
  for(int j = 0; j < n; j++){

    for(int i = 0; i < m; i++){

      int ij = i + m*j;

      std::fill(Tmp.begin(),Tmp.end(),0);

      for(int d = 0; d < ndims; d++){

        Tmp(i,d) = B(j,d) - A(i,d);
      }

      for(int k = 0; k < r; k++){

        X(ij,k) = 0.0;

        for(int ii = 0; ii < m; ii++){

          for(int d = 0; d < ndims; d++){

            int iid = ii + d*m;
            X(ij,k) += Tmp(ii,d)*Q(iid,k);
          }
        }
      }
      
    }
  }

  return X;
}

RcppExport SEXP d_eta_d_phibeta  (SEXP A_, SEXP B_, SEXP Q_){

  using namespace Rcpp;
  NumericMatrix A(A_), B(B_), Q(Q_);

  if (A.ncol() != B.ncol()) return R_NilValue;

  int m = A.nrow();
  int n = B.nrow();
  int ndims = A.ncol();
  int m_ndims = m*ndims;

  if(Q.nrow()!=m_ndims) return R_NilValue;
  int r = Q.ncol();

  NumericMatrix X(n*m,r);
  NumericVector Tmp(m_ndims);

  for(int j = 0; j < n; j++){

    for(int i = 0; i < m; i++){

      int ij = i + m*j;

      std::fill(Tmp.begin(),Tmp.end(),0);
      
      for(int d = 0; d < ndims; d++){

        Tmp(i,d) = B(j,d) - A(i,d);
      }

      for(int k = 0; k < r; k++){

        X(ij,k) = 0.0;

        for(int ii = 0; ii < m; ii++){

          for(int d = 0; d < ndims; d++){

            int iid = ii + d*m;
            X(ij,k) += Tmp(ii,d)*Q(iid,k);
          }
        }
      }

      for(int d = 0; d < ndims; d++){

        X(ij,r + d) = A(i,d);
      }

    }

  }

  return X;
}


