## Practical
rm(list=ls())
set.seed(07062019)
library(survival)
source("~/Dropbox/phd/Package/CaseTimeControl/R/hazardEst.R")
source("~/Dropbox/phd/Package/CaseTimeControl/R/ctc.R")

## Parameters
n <- 1e5
lambda <- rep(131*1e-5,n)
## lambda <- rep(.05/13,n)
tau <- 13
D1 <- 6
D2 <- 3/4
gamma <- log(1)
theta <- log(2)
## lambda <- lambda*exp(c(rep(0,n/2),rep(theta,n/2)))
nSim <- 1e4

## Simulation
estLogit <- or1 <- gammaEst <- numeric(nSim)
for(i in 1:nSim){
    print(i)
    L <- runif(n,0,25)
    ##L <- rexp(n,1/13)
    ## L <- c(runif(n/2,0,25),runif(n/2,0,15))
    u <- runif(n)
    br1 <- 1-exp(-lambda*L)
    br2 <- 1-exp(-lambda*L-lambda*(L+D1)*exp(gamma)+lambda*L*exp(gamma))
    X <- numeric(n)
    X[u<br1] <- -log(1-u[u<br1])/lambda[u<br1]
    X[br1<=u&u<br2] <- (-log(1-u[br1<=u&u<br2])+lambda[br1<=u&u<br2]*L[br1<=u&u<br2]*(exp(gamma)-1))/(lambda[br1<=u&u<br2]*exp(gamma))
    X[br2<=u] <- (-log(1-u[br2<=u])+lambda[br2<=u]*D1*(1-exp(gamma)))/lambda[br2<=u]
    ## X <- rgamma(n, scale=1200,shape=.5)
    ## X <- rlnorm(n,meanlog = log(100))
    cond <- (L<pmin(X,tau)&X>(D1))
    T <- pmin(X,tau)[cond]
    status <- (X<=tau)[cond]
    L <- pmin(L,tau)[cond]
    tmp <- hazardEst(T,status,L,D1,D2,type="short",cut=TRUE)
    ## tmp <- ctc(T,status,L,D1,D2,leftTruncation=TRUE)
    gammaEst[i] <- tmp$Estimate
    or1[i] <- tmp$SE
    print(mean(gammaEst[1:i]))
}
round(exp(mean(gammaEst)),3)
round(mean((gammaEst-qnorm(.975)*or1)<log(1)&(gammaEst+qnorm(.975)*or1)>log(1)),3)

exp(mean(gammaEst) + c(0,-1.96,1.96) * sd(gammaEst) / sqrt(nSim))
mean(or1) + c(0,-1.96,1.96) * sd(or1) / sqrt(nSim)
mean((gammaEst-qnorm(.975)*or1)<log(2)&(gammaEst+qnorm(.975)*or1)>log(1))















