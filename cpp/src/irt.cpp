#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec irt_like(arma::mat X, arma::mat mX, arma::mat pr, arma::vec wts) 
	{
	  arma::vec res = exp(X * log(pr.t()) + mX * log(1-pr.t())) * wts;
	  return res;
	}

