library(survival)

est <- function(time,status,treatmentstart,endofstudy,treatmentlength,window,data=NULL,type="long"){
    if(any(treatmentstart > endofstudy)) warning("Treatment start after end of study")
    time <- eval(substitute(time),data)
    status <- eval(substitute(status),data)
    treatmentstart <- eval(substitute(treatmentstart),data)
    expCases <- (status==1)
    if(!(type%in%c("long","short"))) expCases <- (status==1&time>(treatmentlength+window))
    n <- sum(expCases)
    caseT <- time[expCases]
    caseL <- treatmentstart[expCases]
    casesExpRef <- ((caseT-window)>caseL&(caseT-window)<(caseL+treatmentlength))
    casesExpEvent <- (caseT>caseL&(caseT-treatmentlength)<caseL)
    pBad <- mean(casesExpEvent&!casesExpRef)
    pGood <- mean(!casesExpEvent&casesExpRef)
    LCtrl <- treatmentstart[status==0]
    Fhat <- mean(LCtrl<(endofstudy-treatmentlength))
    if(type=="short") gamma <- log(pBad) - log(pGood) + log(Fhat)
    else if(type=="long"){
        prob1 <- mean(LCtrl>window)
        prob2 <- mean(D2<LCtrl&LCtrl<=(endofstudy-treatmentlength))
        gamma <- log(pBad) - log(pGood) - log(prob1) + log(prob2)
        n2 <- length(LCtrl)
        varZW <- matrix(c(n2*prob2*(1-prob2),n2*(prob2-prob1*prob2),
                          n2*(prob2-prob1*prob2),n2*prob1*(1-prob1)),nrow=2,ncol=2)
    }
    else gamma <- log(pBad)-log(pGood)+log(Fhat)-log(mean(LCtrl>treatmentlength))
    varXY <- matrix(c(n*pBad*(1-pBad),-n*pBad*pGood,-n*pBad*pGood,n*pGood*(1-pGood)),
                    nrow=2,ncol=2)
    nablaf <- c(1/(n*pBad),-1/(n*pGood))
    if(type=="short") Vgamma <- (t(nablaf)%*%varXY%*%nablaf)[1,1] + Fhat*(1-Fhat)/(length(LCtrl)*Fhat^2)
    else if(type=="long"){
        nabla2 <- c(1/(n2*prob2),-1/(n2*prob1))
        Vgamma <- (t(nablaf)%*%varXY%*%nablaf)[1,1] + (t(nabla2)%*%varZW%*%nabla2)[1,1]
    }
    else Vgamma <- (t(nablaf)%*%varXY%*%nablaf)[1,1]
    list(Estimate = gamma, Variance = Vgamma)
}





