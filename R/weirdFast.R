## Practical
rm(list=ls())
set.seed(07062019)
library(survival)
library(foreach)
library(doMC)
registerDoMC(5)
source("~/Dropbox/phd/CTC/R/estimator.R")

## Parameters
n <- 1e5
lambda <- rep(131*1e-5,n)
tau <- 13
D1 <- 6
D2 <- 3/4
gamma <- log(5)
theta <- log(2)
## lambda <- lambda*exp(c(rep(0,n/2),rep(theta,n/2)))
nSim <- 100

## Simulation
estLogit <- or1 <- gammaEst <- numeric(nSim)
sim <- foreach(i = 1:nSim,.combine = "rbind") %dopar% {
    L <- runif(n,0,25)
    ## L <- rexp(n,1/13)
    ## L <- c(runif(n/2,0,25),runif(n/2,0,15))
    u <- runif(n)
    br1 <- 1-exp(-lambda*L)
    br2 <- 1-exp(-lambda*L-lambda*(L+D1)*exp(gamma)+lambda*L*exp(gamma))
    X <- numeric(n)
    X[u<br1] <- -log(1-u[u<br1])/lambda[u<br1]
    X[br1<=u&u<br2] <- (-log(1-u[br1<=u&u<br2])+lambda[br1<=u&u<br2]*L[br1<=u&u<br2]*(exp(gamma)-1))/(lambda[br1<=u&u<br2]*exp(gamma))
    X[br2<=u] <- (-log(1-u[br2<=u])+lambda[br2<=u]*D1*(1-exp(gamma)))/lambda[br2<=u]
    ## X <- rlnorm(n,meanlog = log(100))
    T <- pmin(X,tau)[L<tau]
    status <- (X<=tau)[L<tau]
    L <- L[L<tau]
    tmp <<- est(T,status,L,tau,D1,D2,long=TRUE)
    gammaEst <- tmp$Estimate
    or1 <- tmp$Variance
    c(gammaEst,or1)
}
apply(sim)
exp(mean(or1) + c(0,-1.96,1.96) * sd(or1) / sqrt(nSim))
exp(mean(gammaEst) + c(0,-1.96,1.96) * sd(gammaEst) / sqrt(nSim))
















