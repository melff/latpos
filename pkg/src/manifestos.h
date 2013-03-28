/*#include <Rcpp.h>*/
#include <RcppArmadillo.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
RcppExport SEXP latpos_p (SEXP A_, SEXP B_);
RcppExport SEXP ll_p (SEXP P_, SEXP Y_, SEXP N_, SEXP W_);
RcppExport SEXP d_eta_d_phi (SEXP A_, SEXP B_, SEXP Q_);
RcppExport SEXP d_eta_d_phibeta  (SEXP A_, SEXP B_, SEXP Q_);
RcppExport SEXP latpos_resid (SEXP P_, SEXP Y_, SEXP N_, SEXP W_);
RcppExport SEXP latpos_XWX (SEXP X_, SEXP P_, SEXP N_, SEXP W_);
