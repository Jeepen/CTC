rm(list=ls())
set.seed(31052019)

n <- 1e6
lambda <- 1.31*1e-3
delta1 <- 6
delta2 <- 3/4
C <- 13
C1 <- runif(n,0,25)
C2 <- C1+6
TStar <- rexp(n,lambda)
delta <- (TStar <= C)
T <- pmin(TStar,C)
ExpCases <- (delta==1&C1<T&delta2<TStar)               ## Exposed cases (EC)
nCases <- sum(ExpCases)
casesT <- T[ExpCases]                                  ## Failure time for EC
casesC1 <- C1[ExpCases]                                ## Treatment start for EC
casesC2 <- casesC1 + delta1                            ## End of treatment for EC
## Exposed at reference EC (1/0)
casesExpRef <- ((casesT-delta2)>casesC1&(casesT-delta2)<casesC2)
casesExpEvent <- (casesT>casesC1&casesT<casesC2)      ## Exposed at event EC
controlsExp <- (delta==0&C1<T)                        ## Exposed controls
controlsTStart <- C1[controlsExp][1:nCases]           ## Treatment start for matched exposed controls (MEC)
controlsTEnd <- controlsTStart + delta1               ## End of exposure for MEC
controlsExpRef <- ((casesT-delta2)>controlsTStart&(casesT-delta2)<controlsTEnd)
controlsExpEvent <- (casesT>controlsTStart&casesT<controlsTEnd)
 





