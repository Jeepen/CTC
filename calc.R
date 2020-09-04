## Parameters
tau <- 13
D1 <- 6
D2 <- .75
lambda <- 1e-3
gamma <- log(.5)

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

## Calculations
# Cases -----------------------------------------------------------------------------------------------------------
f1 <- function(t) integrate(function(u) dunif(u,0,25) * lambda * exp(gamma) * exp(-u*lambda - (t-u)*lambda*exp(gamma)), 
                       t-D2, t)$value
f1 <- Vectorize(f1)
badcase <- integrate(f1, 0, tau)$value

f2 <- function(t) integrate(function(u) dunif(u,0,25) * lambda * exp(-D1 * lambda * exp(gamma) - (t-D1) * lambda), 
                            t-D1-D2,t-D1)$value
f2 <- Vectorize(f2)
goodcase <- integrate(f2, 0, tau)$value
(ORcase <- badcase / goodcase)


# Controls --------------------------------------------------------------------------------------------------------
PTL <- function(t){     # P(T=t,L<t)
  if(t > D1) integrate(function(s) lambda * exp(-D1*lambda*exp(gamma)-(t-D1)*lambda) * dunif(s, 0, 25), 0, t-D1)$value + 
    integrate(function(s) lambda*exp(gamma)*exp(-s*lambda-(t-s)*lambda*exp(gamma))*dunif(s,0,25), t-D1, t)$value
  else integrate(function(s) lambda*exp(gamma)*exp(-s*lambda-(t-s)*lambda*exp(gamma))*dunif(s,0,25), 0, t)$value
} 
PTL <- Vectorize(PTL)

PTL2 <- function(t){    # P(T>t,L<t)
  inner <- function(x){ # P(T=x(>t), L<t)
    if(x > D1 & (x-D1)<t) integrate(function(s) lambda*exp(-D1*lambda*exp(gamma)-(x-D1)*lambda)*
                                      dunif(s, 0, 25),0,x-D1)$value + 
      integrate(function(s) lambda*exp(gamma)*exp(-s*lambda-(x-s)*lambda*exp(gamma))*dunif(s,0,25),x-D1, t)$value
    else if(x < D1) integrate(function(s) lambda*exp(gamma)*
                                exp(-s*lambda-(t-s)*lambda*exp(gamma))*dunif(s,0,25), 0, t)$value
    else integrate(function(s) lambda*exp(-D1*lambda*exp(gamma)-(x-D1)*lambda)*
                     dunif(s, 0, 25), 0, t)$value
  }
  inner <- Vectorize(inner)
  integrate(inner, t, Inf)$value
} 

PTL3 <- function(t){    # P(t-D2<L<t,T>t)
  inner <- function(x){ # P(t-D2<L<t,T=x)
    if((x-D1)>(t-D2)&(x-D1)<t) integrate(function(s) lambda*exp(-D1*lambda*exp(gamma)-(x-D1)*lambda)*
                                           dunif(s, 0, 25),t-D2,x-D1)$value + 
      integrate(function(s) lambda*exp(gamma)*exp(-s*lambda-(x-s)*lambda*exp(gamma))*dunif(s,0,25),x-D1,t)$value
    else if((x-D1)>t){
      integrate(function(s) lambda*exp(-D1*lambda*exp(gamma)-(x-D1)*lambda)*
                  dunif(s, 0, 25),t-D2,t)$value
    }
    else{
      integrate(function(s) lambda*exp(gamma)*exp(-s*lambda-(x-s)*lambda*exp(gamma))*dunif(s,0,25),t-D2,t)$value
    }
  }
  inner <- Vectorize(inner)
  integrate(inner, t, Inf)$value
}

PTL4 <- function(t){
  inner <- function(x) lambda*exp(-D1*lambda*exp(gamma)-(x-D1)*lambda) * (punif(t-D1,0,25) - punif(t-D1-D2,0,25))
  inner <- Vectorize(inner)
  integrate(inner, t, Inf)$value
}

badctrl <- integrate(Vectorize(function(t) PTL(t) * PTL3(t) / PTL2(t)), 0, tau)$value
goodctrl <- integrate(Vectorize(function(t) PTL(t) * PTL4(t) / PTL2(t)), 0, tau)$value
(ORctrl <- badctrl / goodctrl)
(OR <- ORcase / ORctrl)
round(log(OR),2)
