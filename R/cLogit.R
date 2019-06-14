## Practical
rm(list=ls())
library(survival)
set.seed(09062019)

## Parameters
n <- 1e6
tau <- 13
D1 <- 6
D2 <- 3/4
lambda <- 131*1e-5

## Function
myLoglik <- function(time,x){
    
}

## Simulation
X <- rexp(n,lambda)
T <- pmin(X,tau)
status <- (X<=tau)
L <- runif(n,0,25)
cases <- (status==1&L<T&D2<T)
casesT <- T[cases]
casesL <- L[cases]
casesExpRef <- (casesT-D1-D2<casesL&casesL<casesT-D1)
casesExpEvent <- (casesT-D2<casesL&casesL<casesT)
nCases <- sum(cases)
dCase <- data.frame(id=rep(1:nCases,2),time=rep(c(1,0),each=nCases),Exp=c(casesExpEvent,casesExpRef))
mCase <- clogit(time~Exp+strata(id),data=dCase)








