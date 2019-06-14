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
n <- 1e5
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
    ExpCases <- (status==1&L<T)#&D2<X)                       ## Exposed cases (EC)
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
    controlsExp <- numeric(nCases)
    for(j in 1:nCases) controlsExp[j] <- sample(which(status==0&L<casesT[j]),1)
    controlsL <- L[controlsExp]
    controlsExpRef <- ((casesT-D2)>controlsL&(casesT-D2)<(controlsL+D1))
    controlsExpEvent <- (casesT>controlsL&(casesT-D1)<controlsL)
    dControl <- data.frame(id=id+nCases,status=rep(c(1,0),each=nCases),
                           Exp=c(controlsExpEvent,controlsExpRef))
    or2[i] <- clogit(status~Exp+strata(id),data=dControl)$coefficients
    print(mean(or2[1:i]))
    ## d <- rbind(dCase,dControl)
    ## d <- transform(d, CaseControl = rep(c(1,0),each=nrow(dCase)))
    ## estLogit[i] <- clogit(status~Exp*CaseControl+strata(id),data=d)$coefficients[3]
    estLogit[i] <- or1[i]-or2[i]
    print(mean(estLogit[1:i]))
}
exp(mean(or1) + c(0,-1.96,1.96) * sd(or1) / sqrt(nSim))
exp(mean(or2) + c(0,-1.96,1.96) * sd(or2) / sqrt(nSim))
exp(mean(estLogit) + c(0,-1.96,1.96)*sd(estLogit)/sqrt(nSim))















