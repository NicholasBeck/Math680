####################################
# Newton Raphson method to solve KERE
# with squared losss in Tweedie model
####################################
setwd("C:/Users/Nicholas/Dropbox/Morvis/Project")
setwd("~/Math680/Code")

library(Rcpp)
library(RcppArmadillo)
library(tweedie)
library(MASS)

source("aobjects.R")
source("kernels.R")
source("kernelmatrix.R")

TweedieQNRKR <- function(X, Y, lambda=1, rho=1.5, phi=1, sig=1){
  ## Setup the kernel matrix
  kern <- rbfdot(sigma=sig)
  Kmat <- kernelMatrix(kern,X) ## symemtric matrix, k(x_i,x_j) = k(x_j,x_i)
  
  ## Setup more values
  n <- nrow(as.matrix(X)) ## Row, Col
  alpha <- rep(1/n,n)
  
  ## CALL THE C CODE to evaluat fhats
  Kmat%*%QNRCpp(Kmat,alpha,Y,rho,lambda,n)
}

sourceCpp('QuasiNewtonRcpp.cpp')

TweedieQNRKR(X,Y)
system.time(TweedieQNRKR(X,Y))
