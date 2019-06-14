## Practical
rm(list=ls())
set.seed(07062019)
library(survival)
source("~/Dropbox/phd/MH.R")

## Parameters
n <- 1e5
lambda <- 131*1e-5
tau <- 13
D1 <- 6
D2 <- 3/4
gamma <- log(2)
theta <- log(2)
lambda <- lambda*exp(c(rep(0,n/2),rep(theta,n/2)))
nSim <- 1000

## Simulation
estLogit <- or1 <- gammaEst <- numeric(nSim)
for(i in 1:nSim){
    print(i)
    ## L <- runif(n,0,25)
    L <- rexp(n,1/13)
    ## L <- c(runif(n/2,0,25),runif(n/2,0,15))
    u <- runif(n)
    br1 <- 1-exp(-lambda*L)
    br2 <- 1-exp(-lambda*L-lambda*(L+D1)*exp(gamma)+lambda*L*exp(gamma))
    ## X <- numeric(n)
    ## X[u<br1] <- -log(1-u[u<br1])/lambda[u<br1]
    ## X[br1<=u&u<br2] <- (-log(1-u[br1<=u&u<br2])+lambda[br1<=u&u<br2]*L[br1<=u&u<br2]*(exp(gamma)-1))/(lambda[br1<=u&u<br2]*exp(gamma))
    ## X[br2<=u] <- (-log(1-u[br2<=u])+lambda[br2<=u]*D1*(1-exp(gamma)))/lambda[br2<=u]
    X <- rlnorm(n,meanlog = log(100))
    T <- pmin(X,tau)
    status <- (X<=tau)
    ExpCases <- (status==1&L<T&(D1+D2)<X)                  ## Exposed cases (EC)
    casesT <- T[ExpCases]                                  ## Failure time for EC
    casesL <- L[ExpCases]                                  ## Treatment start for EC
    casesExpRef <- ((casesT-D2)>casesL&(casesT-D2)<(casesL+D1))
    casesExpEvent <- (casesT>casesL&(casesT-D1)<casesL)    ## Exposed at event EC
    nCases <- length(casesT)                               ## Number of EC
    id <- rep(1:nCases,2)                                  ## Their ID number twice
    ## Collect in data.frame
    dCase <- data.frame(id=id,status=rep(c(1,0),each=nCases),
                        Exp=c(casesExpEvent,casesExpRef))
    or1[i] <- clogit(status~Exp+strata(id),data=dCase)$coefficients
    print(mean(or1[1:i]))
    LCtrl <- L[X>tau]
    gammaEst[i] <- or1[i]-log((mean(LCtrl<tau)-mean(LCtrl<D1))/mean(LCtrl<(tau-D1)))
    ## gammaEst[i] <- or1[i]-log((mean(LCtrl<tau))/mean(LCtrl<(tau-D1)))
    ## gammaEst[i] <- exp(or1)*mean(LCtrl<(tau-D1))/(mean(LCtrl<tau)-mean(LCtrl<D2))
    ## gammaEst[i] <- (or1[i]+D1+log(sum(LCtrl<=(tau-D1))/(sum(LCtrl<=tau)-sum(LCtrl<=D1))))/(1+D1)
    ## gammaEst[i] <- uniroot(function(x) exp(x+D1*(exp(x)-1))*sum(LCtrl<=(tau-D1))/(sum(LCtrl<=tau)-sum(LCtrl<=D1))-exp(or1[i]), interval=c(.01,5))$value
    print(mean(gammaEst[1:i]))
}
exp(mean(or1) + c(0,-1.96,1.96) * sd(or1) / sqrt(nSim))
exp(mean(gammaEst) + c(0,-1.96,1.96) * sd(gammaEst) / sqrt(nSim))

















