## Practical
rm(list=ls())
set.seed(23052019)
library(survival)
library(data.table)
library(foreach)
library(doMC)
library(heaven)

## Parameters
nSim <- 100                                            ## Number of simulations
C <- 13                                                ## 2017-2030
beta <- 0#log(2)                                         ## Assume no association
delta1 <- 6                                            ## Length of treatment
delta2 <- 8/12                                         ## Time between event and reference period
Ntot <- 1e5                                            ## Number of patients
lambda <- rep(1.31*1e-3,Ntot)
## lambda <- c(rep(1.31*1e-3,Ntot/2),rep(2.62*1e-3,Ntot/2))
## gamma <- rep(0,Ntot)                                ## Effect of severity
gamma <- c(rep(0,Ntot/2),rep(log(2),Ntot/2))
lambda <- lambda*exp(gamma)
registerDoMC(5)
estLogit <- estNCC <- numeric(nSim)                    ## Initialize

tStar <- numeric(Ntot)
startTime <- Sys.time()
result <- foreach(i = 1:nSim, .combine = "rbind") %dopar% {
    print(i)
    tStart <- c(runif(Ntot/2,0,25), runif(Ntot/2,0,15))
    tEnd <- tStart + delta1
    startF <- 1 - exp(-lambda*tStart)
    endF <- 1 - exp(-lambda*tStart - (tEnd-tStart)*lambda*exp(beta))
    u <- runif(Ntot)
    cond1 <- u<startF
    cond2 <- u>startF & u<endF
    cond3 <- u>endF
    tStar[cond1] <- -log(1-u[cond1]) / (lambda[cond1])
    tStar[cond2] <- (tStart[cond2]*lambda[cond2]*(exp(beta)-1)-log(1-u[cond2])) /
        (lambda[cond2]*exp(beta))
    tStar[cond3] <- (lambda[cond3]*(tEnd[cond3]-tStart[cond3])*(1-exp(beta))-
                     log(1-u[cond3]))/lambda[cond3]
    T <- pmin(tStar,C)
    delta <- (tStar<=C)                                    ## Status indicator
    ## Cases
    ExpCases <- (delta==1&tStart<T&delta2<tStar)           ## Exposed cases (EC)
    casesT <- T[ExpCases]                                  ## Failure time for EC
    casesTStart <- tStart[ExpCases]                        ## Treatment start for EC
    casesTEnd <- casesTStart + delta1                      ## End of treatment for EC
    ## Exposed at reference EC (1/0)
    casesExpRef <- ((casesT-delta2)>casesTStart&(casesT-delta2)<casesTEnd)
    casesExpEvent <- (casesT>casesTStart&casesT<casesTEnd) ## Exposed at event EC  
    ## casesExpEvent <- (casesExpEvent > casesExpRef)      ## Cheating
    nCases <- length(casesT)                               ## Number of EC
    id <- rep(1:nCases,2)                                  ## Their ID number twice
    ## Collect in data.frame
    dCase <- data.frame(id=id,status=rep(c(1,0),each=nCases),
                        Exp=c(casesExpEvent,casesExpRef))
    ## Conditional logistic regression model for cases
    ## clogit(status~Exp+strata(id),data=dCase)$coefficients
    controlsExp <- numeric(nCases)
    for(j in 1:nCases) controlsExp[j] <- sample(which(T>=casesT[j]&tStart<casesT[j]),1)
    ## controlsExp <- (delta==0&tStart<T)                    ## Exposed controls
    controlsTStart <- tStart[controlsExp][1:nCases]       ## Treatment start for matched exposed controls (MEC)
    controlsTEnd <- controlsTStart + delta1               ## End of exposure for MEC
    ## Exposed at reference
    controlsExpRef <- ((casesT-delta2)>controlsTStart&(casesT-delta2)<controlsTEnd)
    ## Exposed at event
    controlsExpEvent <- (casesT>controlsTStart&casesT<controlsTEnd)
    ## controlsExpEvent <- (controlsExpEvent>controlsExpRef)       ## Cheating
    ## Data.frame for model
    dControl <- data.frame(id=id+nCases,status=rep(c(1,0),each=nCases),
                           Exp=c(controlsExpEvent,controlsExpRef))
    d <- rbind(dCase,dControl)
    d <- transform(d, CaseControl = rep(c(1,0),each=nrow(dCase)))
    ## Conditional logistic regression model for controls
    clogit(status~Exp+strata(id),data=dControl)$coefficients
    ## est[i] <- mCase - mControl                            ## Estimate
    ## clogit(status~Exp*CaseControl+strata(id),data=d)$coefficients[3]
    #----------------------------NCC------------------------------#
    ## nccData <- nccSampling(pnr,time,status,Ncontrols=1,
    ##                        data=data.frame(pnr=1:Ntot,time=T,status=delta))
    ## pnrMatch <- match(nccData$pnr,1:Ntot)
    ## nccData$Exp <- (tStart[pnrMatch] < nccData$time & tEnd[pnrMatch] > nccData$time)
    ## c(clogit(status~Exp*CaseControl+strata(id),data=d)$coefficients[3], clogit(status ~ Exp + strata(time),data=nccData)$coefficients)
}
endTime <- Sys.time()
difftime(endTime,startTime)
apply(result,2,summary)
apply(result,2,function(x)exp(mean(x)+c(0,-1.96,1.96)*sd(x)/sqrt(nSim)))








