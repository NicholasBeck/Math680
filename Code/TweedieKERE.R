####################################
# Newton Raphson method to solve KERE
# with squared losss in Tweedie model
####################################
setwd("~/Math680/Code")
library(Rcpp)
library(RcppArmadillo)
library(tweedie)
library(MASS)

source("aobjects.R")
source("kernels.R")
source("kernelmatrix.R")
sourceCpp('NewtonRcpp.cpp')

TweedieKERE <- function(X, Y, lambda=1, rho=1.5, phi=1, sig=1, n){
  ## Setup the kernel matrix
  kern <- rbfdot(sigma=sig)
  Kmat <- kernelMatrix(kern,X) ## symemtric matrix, k(x_i,x_j) = k(x_j,x_i)
  
  ## Setup more values
  n <- nrow(X) ## Row, Col
  alpha <- seq(0, 1, length.out=n)
  
  ## CALL THE C CODE
  NRCpp(Kmat,alpha,Y,rho,lambda,n)

}
TweedieKERE(X,Y)
system.time(TweedieKERE(X,Y))
