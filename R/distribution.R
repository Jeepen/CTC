## Practical
rm(list=ls())
set.seed(07062019)
library(survival)
source("~/Dropbox/Universitet/phd/MH.R")

## Parameters
lambda <- 131*1e-5
gamma <- log(2)
n <- 1e6
C <- tau <- 13
D1 <- 6
D2 <- 3/4

## Simulation distribution 
X <- rexp(n,lambda)
T <- pmin(X,C)
status <- (X<=C)
L <- runif(n,0,25)
gamma <- log(2)
cases <- (status==1&L<T&D2<T)
ctrl <- (status==0&L<C)

numInt <- function(f,lwr,upr){
    x <- seq(lwr,upr,by=.001)
    sum(f(x)*.001)
}

exp(gamma)*punif(C,0,25)/(punif(C-D1,0,25))

numInt(function(t) numInt(function(s,t) exp(lambda*s*(2-1)),t-D2,t)/25,D2,C)

## Bad cases
(PCase <- integrate(function(t) t/25*dexp(t,lambda),D2,C)$value)
mean(cases)+c(-1.96,1.96)*sd(cases)/sqrt(n)
(PBadCase <- integrate(function(t) dexp(t,lambda)*(punif(t,0,25)-punif(t-D2,0,25)),D2,C)$value)
mean(status==1&L<T&L>T-D2) + c(-1.96,1.96)*sd(status==1&L<T&L>T-D2)/sqrt(n)

## Good cases
(PGoodCase <- integrate(function(t) dexp(t,lambda)*(punif(t-D1,0,25)-punif(t-D1-D2,0,25)),D2,C)$value)
mean(status==1&L<T-D1&L>T-D1-D2) + c(-1.96,1.96)*sd(status==1&L<T-D1&L>T-D1-D2)/sqrt(n)

## Bad ctrl
(PCtrl <- exp(-lambda*C)*C/25)
mean(status==0&L<C)
(PBadCtrl <- integrate(function(t) dexp(t,lambda)*punif(t,0,25)*
                                   (punif(t,0,25)-punif(t-D2,0,25))/punif(C,0,25),D2,C)$value)
(PGoodCtrl <- integrate(function(t) dexp(t,lambda)*punif(t,0,25)*
                                    (punif(t-D1,0,25)-punif(t-D1-D2,0,25))/punif(C,0,25),D2,C)$value)
(PBadCase/PGoodCase)/(PBadCtrl/PGoodCtrl)

## Good ctrl
(PBadCaseBig <- exp(gamma)*integrate(function(t) exp(-t*lambda*exp(gamma))*(punif(t,0,25)-punif(t-D2,0,25)),D2,tau)$value)
(PBadCase <- exp(gamma)*D2*(punif(tau,0,25)-punif(D2,0,25))) 
(PGoodCaseBig <- integrate(function(t) exp(-t*lambda)*(punif(t-D1,0,25)-punif(t-D1-D2,0,25)),D2,tau)$value)
(PGoodCase <- D2*punif(tau-D1,0,25))
PBadCaseBig/PGoodCaseBig
PBadCase/PGoodCase

## OR
(ORcase <- PBadCaseGivenCase/PGoodCaseGivenCase)




ExpCases <- (status==1&L<T&D2<T)           ## Exposed cases (EC)
casesT <- T[ExpCases]                                  ## Failure time for EC
casesTStart <- L[ExpCases]                        ## Treatment start for EC
casesTEnd <- casesTStart + D1                      ## End of treatment for EC
## Exposed at reference EC (1/0)
casesExpRef <- ((casesT-D2)>casesTStart&(casesT-D2)<casesTEnd)
casesExpEvent <- (casesT>casesTStart&casesT<casesTEnd) ## Exposed at event EC
nCases <- length(casesT)                               ## Number of EC
id <- rep(1:nCases,2)                                  ## Their ID number twice
## Collect in data.frame
dCase <- data.frame(id=id,status=rep(c(1,0),each=nCases),
                    Exp=c(casesExpEvent,casesExpRef))
MH(casesT,casesTStart,D1,D2)
est1 <- clogit(status~Exp+strata(id),data=dCase) 
est2 <- sum(casesExpEvent*(1-casesExpRef))/sum((1-casesExpEvent)*casesExpRef)
