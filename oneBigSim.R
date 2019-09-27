## Practical
rm(list=ls())
set.seed(26092019)                                      
library(survival)

## Parameters
n <- 1e6                                                                                          ## Data set size
lambda <- 1e-3                                                                                    ## Hazard
tau <- 13                                                                                         ## End of study
D1 <- 6                                                                                           ## Treatment length
D2 <- .75                                                                                         ## Time between event and reference

## Simulation
L <- runif(n,0,25)                                                                                ## Simulate treatment start                           
X <- rexp(n, lambda)                                                                              ## Simulate failure times
cond <- L<pmin(X,tau)                                                                             ## Truncate data
T <- pmin(X,tau)[cond]                                                                            ## Observed failure time
status <- (X<=tau)[cond]                                                                          ## Case/control
L <- L[cond]                                                                                      ## Observed treatment starts
ExpCases <- (status==1&L<T)                                                                       ## Exposed cases (EC)
casesT <- T[ExpCases]                                                                             ## Failure time for EC
casesL <- L[ExpCases]                                                                             ## Treatment start for EC
nCases <- length(casesT)                                                                          ## Number of cases
casesExpRef <- ((casesT-D2)>casesL&(casesT-D2-D1)<casesL)                                         ## Were the cases exposed at reference?
casesExpEvent <- (casesT>casesL&(casesT-D1)<casesL)                                               ## Were they exposed at event?
## New sampling
controls <- numeric(nCases)                                                                       ## We sample one control for each case
for(j in 1:nCases) controls[j] <- sample(which(T>=casesT[j]&L<casesT[j]),1)                       ## Here the actual sampling is done in the proposed way
controlsT <- T[controls]                                                                          ## Control failure times
controlsL <- L[controls]                                                                          ## Control treatment starts
controlsEE <- (controlsL<casesT&(controlsL+D1)>casesT)                                            ## Controls exposed at event?
controlsER <- (controlsL<(casesT-D2)&controlsL>(casesT-D1-D2))                                    ## Controls exposed at reference?
## Old sampling
controls <- which(status == 0 & L < tau)[1:nCases]                                                ## Classic sampling of controls
controlsT <- T[controls]                                                                          ## Event time for controls
controlsL <- L[controls]                                                                          ## Treatment start for controls
controlsEE0 <- (controlsL<casesT&(controlsL+D1)>casesT)                                            ## Controls exposed at event?
controlsER0 <- (controlsL<(casesT-D2)&controlsL>(casesT-D1-D2))                                    ## Controls exposed at reference

table(casesExpEvent,casesExpRef)
table(controlsEE,controlsER)
table(controlsEE0,controlsER0)
