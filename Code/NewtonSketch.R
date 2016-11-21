source("aobjects.R")
source("kernels.R")
source("kernelmatrix.R")

library("tweedie")
library(MASS)
set.seed(930955)
n <- 1000
d <- 3

# dta <- read.csv("train_set.csv")
# X <- as.matrix(dta[,c(3:5,22:29,31:34)])
# Y <- as.vector(dta[,35])
# n <- nrow(dta)

X <- matrix(nrow = n, ncol = d)
# X[ , 1] <- 1
X[ , 1] <- rnegbin(n, 40, 40)+abs(rnorm(n,0,5)) ## Age
X[ , 2] <- rbinom(1000, 1, 0.6) ## Sex
X[ , 3] <- rgamma(1000, 6, 0.01) ##Income

## Response 

eff <- function(X){
  (-2*(X[, 1] - 30)^2 +  2*X[ , 2] + 30*sqrt(X[, 3]))/500
}
eff2 <- function(X){
  -0.003*(X[ ,1] - 30)^2
}
Y <- rtweedie(n, 1.5, exp(eff2(X)), 1)


rho <- 1.5; lambda <- 15# set these for now
kern <- rbfdot(sigma=1)
Kmat <- kernelMatrix(kern,X[,1]) ## symemtric matrix, k(x_i,x_j) = k(x_j,x_i)

# gradient of the log-likelihood + penalty term
gradient <- function(alpha, Kmat) {
  fhats <- as.vector(alpha %*% Kmat)
  ( exp((2-rho)*fhats) - Y*exp((1-rho)*fhats) + 2*lambda*alpha) %*% Kmat
}

# matrix of second partial derivatives of the log-likelihood + penalty term
hess <- function(alpha, Kmat) {
  fhats <- as.vector(alpha %*% Kmat)
  (Kmat %*% diag(as.vector((2-rho)*exp((2-rho)*fhats) - 
                             (1-rho)*Y*exp((1-rho)*fhats))) %*% Kmat) + 2*lambda*Kmat
}

grad = update = rep(1,n)
alpha = seq(0,1,length.out=n)
## Here is Newton-Raphson!
while(sum(grad^2) > 10^(-5)) {
  grad <- gradient(alpha,Kmat)
  update <- grad %*% solve(hess(alpha,Kmat))
  alpha <- alpha - as.vector(update)
}

fhats <- as.vector(Kmat %*% alpha)
plot(X[,1],(Y-exp(fhats))/(exp(fhats*rho/2))); abline(h=0,lty=2,col='red') ## pearson residuals
plot(X[,3],exp(fhats))
plot(X[,3],eff(X))