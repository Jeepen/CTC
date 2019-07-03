## Practical
rm(list=ls())
set.seed(07062019)
library(survival)
source("~/Dropbox/phd/MH.R")

## Parameters
lambda <- 131*1e-5
tau <- 13
D1 <- 6
D2 <- 3/4
gamma <- log(2)
n <- 1e6
nSim <- 100

## Simulation
estLogit <- or1 <- or2 <- numeric(nSim)
for(i in 1:nSim){
    print(i)
    L <- runif(n,0,25)
    ## L <- rexp(n,1/13)
    u <- runif(n)
    br1 <- 1-exp(-lambda*L)
    br2 <- 1-exp(-lambda*L-lambda*(L+D1)*exp(gamma)+lambda*L*exp(gamma))
    X <- numeric(n)
    X[u<br1] <- -log(1-u[u<br1])/lambda
    X[br1<=u&u<br2] <- (-log(1-u[br1<=u&u<br2])+lambda*L[br1<=u&u<br2]*(exp(gamma)-1))/(lambda*exp(gamma))
    X[br2<=u] <- (-log(1-u[br2<=u])+lambda*D1*(1-exp(gamma)))/lambda
    ## X <- rweibull(n,shape=0.1657955)
    T <- pmin(X,tau)
    status <- (X<=tau)
    ExpCases <- (T>10&T<11)
    casesT <- T[ExpCases]                                  ## Failure time for EC
    casesL <- L[ExpCases]                                  ## Treatment start for EC
    casesExpRef <- (10>casesL&9<casesL)
    casesExpEvent <- (11>casesL&10<casesL)    ## Exposed at event EC
    or1[i] <- log(sum(!casesExpRef&casesExpEvent)) - log(sum(casesExpRef&!casesExpEvent))
    print(mean(or1[1:i]))
}
exp(mean(or1) + c(0,-1.96,1.96) * sd(or1) / sqrt(nSim))
exp(mean(or2) + c(0,-1.96,1.96) * sd(or2) / sqrt(nSim))
exp(mean(estLogit) + c(0,-1.96,1.96)*sd(estLogit)/sqrt(nSim))















