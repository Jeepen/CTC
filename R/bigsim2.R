## Practical
rm(list=ls())
set.seed(07062019)
source("~/Dropbox/phd/CTC/R/estimator.R")

## Parameters
n <- 1e6
lambda <- c(1e-2,1e-3,1e-4)
tau <- 13
D1 <- 6
D2 <- 3/4
gamma <- c(log(1),log(2),log(3))
Ls <- c("gamma","lnorm","wei")
nSim <- 100

## Simulation
means <- estSD <- empSD <- matrix(NA,nrow=9,ncol=3)
SE <- gammaEst <- numeric(nSim)
X <- numeric(n)
tab <- matrix(NA, nrow = 3, ncol = 3)
for(x in 1:3){
  for(y in 1:3){
    ## for(z in 1:3){
      ## print((x-1)*9+(y-1)*3+z)
      print(y+(x-1)*3)
      for(i in 1:nSim){
        ## L <- switch(Ls[z], gamma = rgamma(n,shape=10), lnorm = rlnorm(n,meanlog = log(10)-.5),
          ## wei = rweibull(n,shape=1/3.4))
          L <- runif(n,0,20)
        u <- runif(n)
        br1 <- 1-exp(-lambda[x]*L)
        br2 <- 1-exp(-lambda[x]*L-lambda[x]*(L+D1)*exp(gamma[y])+lambda[x]*L*exp(gamma[y]))
        X[u<br1] <- -log(1-u[u<br1])/lambda[x]
        X[br1<=u&u<br2] <- (-log(1-u[br1<=u&u<br2])+lambda[x]*L[br1<=u&u<br2]*(exp(gamma[y])-1))/(lambda[x]*exp(gamma[y]))
        X[br2<=u] <- (-log(1-u[br2<=u])+lambda[x]*D1*(1-exp(gamma[y])))/lambda[x]
        cond <- L<X&X<tau&X>(D1+D2)
        T <- X[cond]      
        L <- L[cond]
        gammaEst[i] <- log(sum(((T-D2)<L&L<T))) - log(sum((T-D1-D2)<L&L<(T-D1)))
      }
      tab[x,y] <- exp(mean(gammaEst))
      print(tab)
      ## means[(x-1)*3+z,y] <- exp(mean(gammaEst))
      ## print(means[(x-1)*3+z,y])
    ## }
  }
}








