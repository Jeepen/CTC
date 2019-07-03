rm(list=ls())

## Functions
PBadCase <- function(D2,tau,lambda) integrate(function(t) dexp(t,lambda)*(punif(t,0,25)-punif(t-D2,0,25)),0,tau)$value
PGoodCase <- function(D1,D2,tau,lambda) integrate(function(t) dexp(t,lambda)*(punif(t-D1,0,25)-punif(t-D1-D2,0,25)), 0, tau)$value
PBadCtrl <- function(D2,tau,lambda) integrate(function(t) dexp(t,lambda)*punif(t,0,25)*(punif(t,0,25)-punif(t-D2,0,25))/punif(tau,0,25), 0, tau)$value
PGoodCtrl <- function(D1,D2,tau,lambda) integrate(function(t) dexp(t,lambda)*punif(t,0,25)*(punif(t-D1,0,25)-punif(t-D1-D2,0,25))/punif(tau,0,25), 0, tau)$value

## Calculations
tau <- 13
lambda <- 1.31*1e-3
D1 <- 6
D2 <- 3/4
(ORCase <- PBadCase(Delta2,tau,lambda)/PGoodCase(Delta1,Delta2,tau,lambda))
(ORCtrl <- PBadCtrl(Delta2,tau,lambda)/PGoodCtrl(Delta1,Delta2,tau,lambda))
ORCase/ORCtrl
