#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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
arma::vec QNRCpp(arma::mat Kmat, arma::vec alpha_start, arma::vec Y, double rho, double lambda, int n){
  //Identity used everywhere
  arma::mat Iden=eye(n,n);
  
  //Base declarations
  arma::vec update(n);
  arma::vec alpha=alpha_start;
  arma::mat dlam = lambda*Iden;
  arma::vec fhats(n);
  arma::mat hess(n,n); 
  
  //Quasi Newton-Pieces
  arma::vec y(n);
  arma::vec s(n);
  double ys;
  arma::vec gradnow(n);
  arma::vec gradold(n);
  arma::mat hessI(n,n);
  
  //First Naive Update
  fhats   = Kmat*alpha;
  gradnow = gradCpp(Kmat,fhats,Y,rho,lambda);
  hess    = hessCpp(Kmat,fhats,Y,rho,dlam);
  hessI   = inv(hess);
  update  = hessI*gradnow;
  alpha  -= update;
  
  //Second Updated (get better starting Hess Inverse)
  fhats   = Kmat*alpha;
  gradnow = gradCpp(Kmat,fhats,Y,rho,lambda);
  hess    = hessCpp(Kmat,fhats,Y,rho,dlam);
  hessI   = inv(hess);
  update  = hessI*gradnow;
  alpha  -= update;
  
  //Resume looping
  double check = 1;
  int i=3;
  while(check > pow(10,-5)) {
  //for(int i=2;i<10;i++){
    gradold = gradnow;
    fhats = Kmat*alpha;
    gradnow = gradCpp(Kmat,fhats,Y,rho,lambda);
    
    s = -update;
    y = gradnow-gradold;
    ys = dot(y,s);
    
    hessI = (Iden-s*y.t()/ys)*hessI*(Iden-y*s.t()/ys)+s*s.t()/ys;
    update=hessI*gradnow;
    alpha -= update;
    
    check = dot(gradnow,gradnow);
    i+=1;
    
    // hessI = (Iden-s*y.t()/ys)*hessI*(Iden-y*s.t()/ys)+s*s.t()/ys;
    // update = hessI*gradnow;
    // alpha -= update;
    // s = -update;
    // y = gradnow - gradold ;
    // ys = dot(y,s);
    // hessI = (Iden-s*y.t()/ys)*hessI*(Iden-y*s.t()/ys)+s*s.t()/ys;
    // gradold = gradnow;
    // 
    // fhats = Kmat*alpha;
    // gradnow = gradCpp(Kmat,fhats,Y,rho,lambda);
    // 
    // check = dot(gradold,gradold);
    // i+=1;
  }
  return alpha;
}
