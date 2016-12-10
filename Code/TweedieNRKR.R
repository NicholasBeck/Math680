####################################
# Newton Raphson method to solve KERE
# with squared losss in Tweedie model
####################################
setwd("C:/Users/Nicholas/Dropbox/Morvis/Project")
library(Rcpp)
library(RcppArmadillo)
library(tweedie)
library(MASS)

source("aobjects.R")
source("kernels.R")
source("kernelmatrix.R")

TweedieNRKR <- function(X, Y, lambda=1, rho=1.5, phi=1, sig=1){
  ## Setup the kernel matrix
  kern <- rbfdot(sigma=sig)
  Kmat <- kernelMatrix(kern,X) ## symemtric matrix, k(x_i,x_j) = k(x_j,x_i)
  
  ## Setup more values
  n <- nrow(X) ## Row, Col
  alpha <- rep(1/n,n)
  
  ## CALL THE C CODE
  Kmat%*%NRCpp(Kmat,alpha,Y,rho,lambda,n)

}

sourceCpp('NewtonRcpp.cpp')

TweedieKERE(X,Y)
system.time(TweedieKERE(X,Y))
