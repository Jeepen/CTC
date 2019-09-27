## Practical
rm(list=ls())
set.seed(26092019)                                      
library(survival)

## Parameters
n <- 1e5                                                                                          ## Data set size
lambda <- 1e-3                                                                                    ## Hazard
tau <- 13                                                                                         ## End of study
D1 <- 6                                                                                           ## Treatment length
D2 <- .75                                                                                         ## Time between event and reference
nsim <- 1000                                                                                      ## Number of simulations

## Simulation
estNew <- timetrend <- estOld <- estID <- numeric(nsim)
for(i in 1:nsim){
  print(i)
  L <- runif(n,0,25)                                                                              ## Simulate treatment start                           
  X <- rexp(n, lambda)                                                                            ## Simulate failure times
  cond <- L<pmin(X,tau)                                                                           ## Truncate data
  T <- pmin(X,tau)[cond]                                                                          ## Observed failure time
  status <- (X<=tau)[cond]                                                                        ## Case/control
  L <- L[cond]                                                                                    ## Observed treatment starts
  ExpCases <- (status==1&L<T)                                                                     ## Exposed cases (EC)
  casesT <- T[ExpCases]                                                                           ## Failure time for EC
  casesL <- L[ExpCases]                                                                           ## Treatment start for EC
  nCases <- length(casesT)                                                                        ## Number of cases
  casesExpRef <- ((casesT-D2)>casesL&(casesT-D2-D1)<casesL)                                       ## Were the cases exposed at reference?
  casesExpEvent <- (casesT>casesL&(casesT-D1)<casesL)                                             ## Were they exposed at event?
  dcase <- data.frame(time = rep(c(0,1),each = nCases),
                      exp = c(casesExpRef, casesExpEvent),
                      id = rep(1:nCases, 2), case = 1)                                            ## Case part of data set for conditional logistic regression
  ## New sampling
  controls <- numeric(nCases)                                                                     ## We sample one control for each case
  for(j in 1:nCases) controls[j] <- sample(which(T>=casesT[j]&L<casesT[j]),1)                     ## Here the actual sampling is done in the proposed way
  controlsT <- T[controls]                                                                        ## Control failure times
  controlsL <- L[controls]                                                                        ## Control treatment starts
  controlsEE <- (controlsL<casesT&(controlsL+D1)>casesT)                                          ## Controls exposed at event?
  controlsER <- (controlsL<(casesT-D2)&controlsL>(casesT-D1-D2))                                  ## Controls exposed at reference?
  dcont <- data.frame(time = rep(c(0,1), each = nCases),
                      exp = c(controlsER, controlsEE),
                      id = rep(nCases+1:nCases, 2), case = 0)                                     ## Control part of data set for conditional logistic regression
  dNew <- rbind(dcase,dcont)                                                                      ## Merge case and control part of data set into one data set
  m <- clogit(time ~ exp + exp*case + strata(id), data = dNew)                                    ## Fit conditional logistic regression model
  estNew[i] <- coef(m)[3]                                                                         ## Return desired estimate
  timetrend[i] <- coef(m)[1]                                                                      ## Odds ratio for cases
  print(mean(estNew[1:i]))
  print(mean(timetrend[1:i]))
  ## Old sampling
  controls <- which(status == 0 & L < tau)[1:nCases]                                              ## Classic sampling of controls
  controlsT <- T[controls]                                                                        ## Event time for controls
  controlsL <- L[controls]                                                                        ## Treatment start for controls
  controlsEE <- (controlsL<casesT&(controlsL+D1)>casesT)                                          ## Controls exposed at event?
  controlsER <- (controlsL<(casesT-D2)&controlsL>(casesT-D1-D2))                                  ## Controls exposed at reference
  dcont <- data.frame(time = rep(c(0,1), each = nCases),
                      exp = c(controlsER, controlsEE),
                      id = rep(nCases+1:nCases, 2), case = 0)                                     ## Control part of data set for conditional logistic regression
  dOld <- rbind(dcase,dcont)                                                                      ## Merge case and control part of data set into one data set
  m <- clogit(time ~ exp + exp*case + strata(id), data = dOld)                                    ## Fit conditional logistic regression
  estOld[i] <- coef(m)[3]                                                                         ## Return desired estimate
  print(mean(estOld[1:i]))
  ## Old sampling wrong ID
  dcont <- data.frame(time = rep(c(0,1), each = nCases),
                      exp = c(controlsER, controlsEE),
                      id = rep(1:nCases, 2), case = 0)                                            ## Control part of data set with wrong id
  dOld <- rbind(dcase,dcont)                                                                      ## Merge case and control part of data set into one data set
  m <- clogit(time ~ exp + exp*case + strata(id), data = dOld)                                    ## Fit conditional logistic regression model
  estID[i] <- coef(m)[3]                                                                          ## Return desired estimate
  print(mean(estID[1:i])) 
}
round(exp(apply(cbind(estNew,timetrend,estOld,estID),2,mean)),2)                                  ## Get odds ratios
## 1.00, 1.92, 1.43, 1.26


