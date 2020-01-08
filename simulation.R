## Seed and packages
set.seed(16122019)
library(survival)
library(sandwich)
library(foreach)
library(mets)
library(doMC)
registerDoMC(6)
getDoParWorkers()

## Parameters
n <- 3e4
lambda <- 1e-2*rgamma(n,1,1)
tau <- 13
D1 <- 6
D2 <- 3/4
gamma <- c(log(1),log(.5),log(2),log(5))
nsim <- 1000

## Simulation
proposed <- classic <- sdEst <- sandwichEst <- rep(NA,nrow=nsim)
results <- foreach(w=1:4) %dopar% {
    for(i in 1:nsim){
        if((i%%10)==0){
            print(w)
            print(i)
        }
        L <- rweibull(n,shape=2,scale=14)
        u <- runif(n)
        br1 <- 1-exp(-lambda*L)
        br2 <- 1-exp(-lambda*L-lambda*(L+D1)*exp(gamma[w])+lambda*L*exp(gamma[w]))
        X <- numeric(n)
        X[u<br1] <- -log(1-u[u<br1])/lambda[u<br1]
        X[br1<=u&u<br2] <- (-log(1-u[br1<=u&u<br2])+
                            lambda[br1<=u&u<br2]*L[br1<=u&u<br2]*(exp(gamma[w])-1))/(lambda[br1<=u&u<br2]*exp(gamma[w]))
        X[br2<=u] <- (-log(1-u[br2<=u])+lambda[br2<=u]*D1*(1-exp(gamma[w])))/lambda[br2<=u]
        T <- pmin(X,tau)
        status <- (X<=tau)
        ExpCases <- (status==1&L<T)                            ## Exposed cases (EC)
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
        controlsExp <- which(status==0,L<tau)[1:nCases]
        controlsL <- L[controlsExp]
        controlsExpRef <- ((casesT-D2)>controlsL&(casesT-D2)<(controlsL+D1))
        controlsExpEvent <- (casesT>controlsL&(casesT-D1)<controlsL)
        dControl <- data.frame(id=id+nCases,clusters=id,status=rep(c(1,0),each=nCases),
                               Exp=c(controlsExpEvent,controlsExpRef), case = 0)
        d <- rbind(dCase,dControl)
        d <- transform(d, EG = Exp*case)
        tmp <- phreg(Surv(rep(1,nrow(d)),status) ~ Exp + EG + strata(id), data=d)
        classic[i] <- coef(tmp)[2]
        ## Sample controls in proposed way
        for(j in 1:nCases) controlsExp[j] <- sample(which(status==0&L<casesT[j]),1)
        controlsL <- L[controlsExp]
        controlsExpRef <- ((casesT-D2)>controlsL&(casesT-D2)<(controlsL+D1))
        controlsExpEvent <- (casesT>controlsL&(casesT-D1)<controlsL)
        dControl <- data.frame(id=id+nCases,clusters=id,status=rep(c(1,0),each=nCases),
                               Exp=c(controlsExpEvent,controlsExpRef), case = 0)
        d <- rbind(dCase,dControl)
        d <- transform(d, EG = Exp*case)
        tmp <- phreg(Surv(rep(1,nrow(d)), status) ~ Exp + EG + strata(id), data=d)
        proposed[i] <- coef(tmp)[2]
        sdEst[i] <- sqrt(vcov(tmp)[2,2])
        ## sdEst[i] <- summary(tmp)$coefficients[2,3]
        tmp <- phreg(Surv(rep(1,nrow(d)), status) ~ Exp + EG + strata(id) + cluster(clusters), data=d)
        sandwichEst[i] <- sqrt(vcov(tmp)[2,2])
        ## Sanwich estimator
        ## helper <- estfun(tmp)
        ## helper <- helper[1:nCases,] + helper[(nCases+1):(2*nCases),] + 
	## helper[(2*nCases+1):(3*nCasesases),] + helper[(3*nCases+1):(4*nCases),]
        ## fyld <- matrix(0,nrow=2,ncol=2)
        ## for(j in 1:nCases) fyld <- fyld + helper[j,]%*%t(helper[j,])
        ## brot <- vcov(tmp)
        ## sandwichEst[i] <- sqrt((brot%*%fyld%*%brot)[2,2])
    }
    data.frame(classic,proposed,sdEst,sandwichEst)
}
result <- do.call("rbind", lapply(1:4, function(x) apply(results[[x]], 2, mean)))
result <- cbind(result, emp = do.call("c", lapply(1:4, function(x) sd(results[[x]][,2]))))


