#include <RcppArmadillo.h>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]


double LOCO_TS(arma::vec union_lambda, arma::mat coefs, arma::mat coefs_j){
  
  double TS; 
  
  int M = union_lambda.n_rows;
  
  arma::mat Delta = coefs_j - coefs;
  arma::mat Delta_1 = Delta.rows(0, M-2);
  arma::mat Delta_2 = Delta.rows(1, M-1);
  arma::vec Lambda = diff(union_lambda);
  
  arma::mat Delta_sq = 1.0/3 * (Delta_1 % Delta_1 + Delta_1 % Delta_2 + Delta_2 % Delta_2 );
  
  arma::mat Epsilon = Delta_sq.each_col() % Lambda;
  
  TS = accu(Epsilon);
  
  
  
  
  return TS;
  
}



