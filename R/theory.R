rm(list=ls())
library(pracma)

## Functions
PBadCase <- function(Delta2,tau,lambda){
    Delta2*integral(function(t) lambda*exp(-lambda*t),Delta2,tau)
}
PGoodCase <- function(Delta1,Delta2,tau,lambda){
    Delta2*integral(function(t) lambda*exp(-lambda*t),Delta1+Delta2,tau)+
        integral(function(t) (t-Delta1)*lambda*exp(-lambda*t),Delta1,Delta1+Delta2)
}
PBadCtrl <- function(Delta2,tau,lambda){
    Delta2*integral(function(t) t*lambda*exp(-lambda*t),Delta2,tau)
}
PGoodCtrl <- function(Delta1,Delta2,tau,lambda){
    Delta2*integral(function(t) t*lambda*exp(-lambda*t),Delta1+Delta2,tau) +
        integral(function(t) t*(t-Delta1)*lambda*exp(-lambda*t),Delta1,Delta1+Delta2)
}

## Calculations
tau <- 13
lambda <- 1.31*1e-3
Delta1 <- 6
Delta2 <- 2/3
(ORCase <- PBadCase(Delta2,tau,lambda)/PGoodCase(Delta1,Delta2,tau,lambda))
(ORCtrl <- PBadCtrl(Delta2,tau,lambda)/PGoodCtrl(Delta1,Delta2,tau,lambda))
ORCase/ORCtrl
