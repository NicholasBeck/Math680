#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec gradCpp(arma::mat Kmat, arma::vec alpha, arma::vec Y, double rho, double lambda){
  arma::vec fhats = Kmat * alpha;
  return Kmat*(exp((2-rho)*fhats)-Y%exp((1-rho)*fhats)+2*lambda*alpha);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
// matrix of second partial derivatives of the log-likelihood + penalty term
arma::mat hessCpp(arma::mat Kmat, arma::vec alpha, arma::vec Y, double rho, double lambda){
  arma::vec fhats = Kmat * alpha;
  return Kmat * diagmat((2-rho)*exp((2-rho)*fhats)-(1-rho)*Y%exp((1-rho)*fhats))* Kmat + 2*lambda*Kmat;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec NRCpp(arma::mat Kmat, arma::vec alpha_start, arma::vec Y, double rho, double lambda){
int n=Kmat.n_rows;
arma::vec grad=arma::ones<arma::vec>(n);
arma::vec update(n);
arma::vec alpha=alpha_start;
double check = 1;
int i=1;
while(check > pow(10,-5)) {
    grad = gradCpp(Kmat,alpha,Y,rho,lambda);
    update = inv(hessCpp(Kmat,alpha,Y,rho,lambda))*grad;
    alpha -= update;
    check = accu(grad%grad);
    i+=1;
  }
return alpha;
}
