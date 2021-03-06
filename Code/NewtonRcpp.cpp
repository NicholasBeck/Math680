#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::vec fh(arma::mat Kmat, arma::vec alpha){
//   return Kmat * alpha;
// }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec gradCpp(arma::mat Kmat, arma::vec fhats, arma::vec Y, double rho, double lambda){
  return Kmat*(exp((2-rho)*fhats)-Y%exp((1-rho)*fhats))+lambda*fhats;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
// matrix of second partial derivatives of the log-likelihood + penalty term
arma::mat hessCpp(arma::mat Kmat, arma::vec fhats, arma::vec Y, double rho, arma::mat dlam){
  return (Kmat * diagmat((2-rho)*exp((2-rho)*fhats)-(1-rho)*Y%exp((1-rho)*fhats)) + dlam)*Kmat;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec NRCpp(arma::mat Kmat, arma::vec alpha_start, arma::vec Y, double rho, double lambda,int n){
arma::vec grad(n);
arma::vec update(n);
arma::vec alpha=alpha_start;
arma::vec fhats(n);
arma::mat dlam = lambda*eye(n,n);
double check = 1;
int i=1;
while(check > pow(10,-5)) {
    fhats = Kmat*alpha;
    grad = gradCpp(Kmat,fhats,Y,rho,lambda);
    update = inv(hessCpp(Kmat,fhats,Y,rho,dlam))*grad;
    alpha -= update;
    check = dot(grad,grad);
    i+=1;
  }
return alpha;
}
