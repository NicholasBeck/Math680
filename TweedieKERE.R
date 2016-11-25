####################################
# Newton Raphson method to solve KERE
# with squared losss in Tweedie model
####################################

library(Rcpp)
source("aobjects.R")
source("kernels.R")
source("kernelmatrix.R")


TweedieKERE <- function(X, Y, lambda=1, rho=1.5, phi=1, sig=1){
  ## Setup the kernel matrix
  kern <- rbfdot(sigma=sig)
  Kmat <- kernelMatrix(kern,X) ## symemtric matrix, k(x_i,x_j) = k(x_j,x_i)
  
  ## Setup more values
  dimsX <- dim(X) ## Row, Col
  alpha <- seq(0, 1, length.out=dimsX[1])
  ## CALL THE C CODE
  
  ##dyn.load('tempName.so') ## 
  ##.Call('tempName', X, Y, Kmat, alpha, lambda, rho, phi, dimsX[1], dimsX[2])

}