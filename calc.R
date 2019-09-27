## Parameters
tau <- 13
D1 <- 6
D2 <- .75
lambda <- 1e-3

## Calculations
(ee <- integrate(function(t) dexp(t,lambda)*(punif(t,0,25)-punif(t-D2,0,25)),0,tau)$value)
(er <- integrate(function(t) dexp(t,lambda)*(punif(t-D1,0,25)-punif(t-D2-D1,0,25)),0,tau)$value)
round(ee/er,2)
