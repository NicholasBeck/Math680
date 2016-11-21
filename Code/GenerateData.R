## Generate training data
install.package("tweedie")
library("tweedie")
library(MASS)
set.seed(930955)
n <- 1000
d <- 3

X <- matrix(nrow = n, ncol = d)
X[ , 1] <- rnegbin(n, 40, 40) ## Age 
X[ , 2] <- rbinom(n, 1, 0.6) ## Sex 
X[ , 3] <- rgamma(n, 6, 0.01) ##Income

## Response 

eff <- function(X){
  (-2*(X[, 1] - 30)^2 +  2*X[ , 2] + 30*sqrt(X[, 3]))/500
}
eff2 <- function(X){
  -0.003*(X[ ,1] - 30)^2
}

Y <- rtweedie(n, 1.5, exp(eff2(X)), 1)

### MOre crap

X <- matrix(nrow = 1000, ncol = 1)
X[, 1] <- rgamma(1000, 6, 0.1)

Y <- rtweedie(1000, 1.5, exp(eff2(X)), 1)

