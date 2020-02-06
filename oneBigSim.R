## Seed and packages
rm(list=ls())
set.seed(16122019)
# set.seed(26092019)

## Parameters
n <- 1e6
lambda <- 1e-3 
tau <- 13
D1 <- 6
D2 <- 3/4

## Simulation
L <- runif(n,0,25)
X <- rexp(n, lambda)    
T <- pmin(X,tau)
ExpCases <- (X<=tau & L<T)                             ## Exposed cases (EC)
casesT <- T[ExpCases]                                  ## Failure time for EC
casesL <- L[ExpCases]                                  ## Treatment start for EC
casesExpRef <- ((casesT-D2)>casesL&(casesT-D2)<(casesL+D1))
casesExpEvent <- (casesT>casesL&(casesT-D1)<casesL)    ## Exposed at event EC
nCases <- length(casesT)                               ## Number of EC
id <- rep(1:nCases,2)                                  ## Their ID number twice
## Collect in data.frame
dCase <- data.frame(id=id,clusters=id,status=rep(c(1,0),each=nCases),
                    Exp=c(casesExpEvent,casesExpRef), case = 1)
## Sample controls in usual way
controlsExp <- which(X>tau,L<tau)[1:nCases]
controlsL <- L[controlsExp]
controlsExpRef <- ((casesT-D2)>controlsL&(casesT-D2)<(controlsL+D1))
controlsExpEvent <- (casesT>controlsL&(casesT-D1)<controlsL)
dControl <- data.frame(id=id+nCases,clusters=id,status=rep(c(1,0),each=nCases),
                       Exp=c(controlsExpEvent,controlsExpRef), case = 0)
## Sample controls in proposed way
for(j in 1:nCases) controlsExp[j] <- sample(which(T>=casesT[j]&L<casesT[j]),1)
controlsL <- L[controlsExp]
controlsExpRef2 <- ((casesT-D2)>controlsL&(casesT-D2)<(controlsL+D1))
controlsExpEvent2 <- (casesT>controlsL&(casesT-D1)<controlsL)

table(casesExpEvent,casesExpRef)
table(controlsExpEvent,controlsExpRef)
table(controlsExpEvent2,controlsExpRef2)

