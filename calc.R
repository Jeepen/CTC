## Parameters
tau <- 13
D1 <- 6
D2 <- .75
lambda <- 1e-3

## Odds ratio for cases
(ee <- integrate(function(t) dexp(t,lambda)*(punif(t,0,25)-punif(t-D2,0,25)),0,tau)$value)
(er <- integrate(function(t) dexp(t,lambda)*(punif(t-D1,0,25)-punif(t-D2-D1,0,25)),0,tau)$value)
(or1 <- ee/er)
round(or1,2)

## Odds ratio for controls when sampling incorrectly
(eec <- integrate(function(t) dexp(t,lambda)*punif(t,0,25)/punif(tau,0,25)*(punif(t,0,25)-punif(t-D2,0,25)),0,tau)$value)
(erc <- integrate(function(t) dexp(t,lambda)*punif(t,0,25)/punif(tau,0,25)*(punif(t-D1,0,25)-punif(t-D2-D1,0,25)),0,tau)$value)
(or2 <- or1 / (eec / erc))
round(or2,2)

