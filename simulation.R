## Seed and packages
rm(list=ls())
library(survival)
library(foreach)
library(mets)
library(doRNG)
library(doMC)
registerDoMC(6)
getDoParWorkers()

## Parameters
n <- 1e5
lambda <- rep(1e-3,n)
tau <- 13
D1 <- 6
D2 <- 3/4
gamma <- c(log(1), log(.5), log(2), log(5))
nsim <- 10000

## Simulation
proposed <- classic <- sdEst <- sandwichEst <- rep(NA,nrow=nsim)
results2 <- foreach(w=1:4, .packages = "mets", .options.RNG = 16122019) %dorng% {
    X <- numeric(n)
    for(i in 1:nsim){
        if((i%%10)==0){
            print(w)
            print(i)
        }
        L <- runif(n, 0, 25)
        u <- runif(n)
        br1 <- 1 - exp(-lambda * L)
        br2 <- 1 - exp(-lambda * L - lambda * (L + D1) * exp(gamma[w]) + lambda * L * exp(gamma[w]))
        X[u < br1] <- -log(1 - u[u < br1]) / lambda[u < br1]
        X[br1 <= u & u < br2] <- (-log(1 - u[br1 <= u & u < br2]) +
                                  lambda[br1 <= u & u < br2] * L[br1 <= u & u < br2] * (exp(gamma[w]) - 1)) /
            (lambda[br1 <= u & u < br2] * exp(gamma[w]))
        X[br2<=u] <- (-log(1 - u[br2 <= u]) + lambda[br2 <= u] * D1 * (1 - exp(gamma[w]))) / lambda[br2 <= u]
        T <- pmin(X, tau)
        ExpCases <- (X <= tau & L < T)                         ## Exposed cases (EC)
        casesT <- T[ExpCases]                                  ## Failure time for EC
        casesL <- L[ExpCases]                                  ## Treatment start for EC
        casesExpRef <- ((casesT - D2) > casesL & (casesT - D2) < (casesL + D1))
        casesExpEvent <- (casesT > casesL & (casesT - D1) < casesL)    ## Exposed at event EC
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
        d <- rbind(dCase,dControl)
        d <- transform(d, EG = Exp*case)
        tmp <- coxph(Surv(rep(1,nrow(d)), status) ~ Exp + EG + strata(id), data=d)
        classic[i] <- coef(tmp)[2]
        ## Sample controls in proposed way
        for(j in 1:nCases) controlsExp[j] <- sample(which(T>=casesT[j]&L<casesT[j]),1)
        controlsL <- L[controlsExp]
        controlsExpRef <- ((casesT-D2)>controlsL&(casesT-D2)<(controlsL+D1))
        controlsExpEvent <- (casesT>controlsL&(casesT-D1)<controlsL)
        dControl <- data.frame(id=id+nCases,clusters=id,status=rep(c(1,0),each=nCases),
                               Exp=c(controlsExpEvent,controlsExpRef), case = 0)
        d <- rbind(dCase,dControl)
        d <- transform(d, EG = Exp*case)
        tmp <- coxph(Surv(rep(1,nrow(d)), status) ~ Exp + EG + strata(id), data=d)
        proposed[i] <- coef(tmp)[2]
        sdEst[i] <- sqrt(vcov(tmp)[2,2])
        tmp <- phreg(Surv(rep(1,nrow(d)), status) ~ Exp + EG + strata(id) + cluster(clusters), data=d)
        sandwichEst[i] <- sqrt(vcov(tmp)[2,2])
    }
    data.frame(classic,proposed,sdEst,sandwichEst)
}
result2 <- do.call("rbind", lapply(1:4, function(x) apply(results2[[x]], 2, mean)))
result2 <- cbind(result2, emp = do.call("c", lapply(1:4, function(x) sd(results2[[x]][,2]))))
result2
round(result2, 3)

mean(results2[[1]][,2]-qnorm(.975)*results2[[1]][,3]<0&results2[[1]][,2]+qnorm(.975)*results2[[1]][,3]>0) # 0.9603
mean(results2[[1]][,2]-qnorm(.975)*results2[[1]][,4]<0&results2[[1]][,2]+qnorm(.975)*results2[[1]][,4]>0) # 0.9525

saveRDS(results2, "~/Dropbox/phd/articles/CTCcode/simResults.rds")
