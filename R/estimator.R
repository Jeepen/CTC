est <- function(time,status,treatmentstart,endofstudy,treatmentlength,window,data=NULL,type="short"){
    if(any(treatmentstart > endofstudy)) warning("Treatment start after end of study")
    time <- eval(substitute(time),data)
    status <- eval(substitute(status),data)
    treatmentstart <- eval(substitute(treatmentstart),data)
    expCases <- (status==1)
    n <- sum(expCases)
    caseT <- time[expCases]
    caseL <- treatmentstart[expCases]
    casesExpRef <- ((caseT-window)>caseL&(caseT-window)<(caseL+treatmentlength))
    casesExpEvent <- (caseT>caseL&(caseT-treatmentlength)<caseL)
    pBad <- mean(casesExpEvent&!casesExpRef)
    pGood <- mean(!casesExpEvent&casesExpRef)
    LCtrl <- treatmentstart[status==0]
    Fhat <- mean(LCtrl<(endofstudy-treatmentlength))
    gamma <- log(pBad) - log(pGood) + log(Fhat)
    varXY <- matrix(c(n*pBad*(1-pBad),-n*pBad*pGood,-n*pBad*pGood,n*pGood*(1-pGood)),
                    nrow=2,ncol=2)
    nablaf <- c(1/(n*pBad),-1/(n*pGood))
    Vgamma <- (t(nablaf)%*%varXY%*%nablaf)[1,1] + Fhat*(1-Fhat)/(length(LCtrl)*Fhat^2)
    list(Estimate = gamma, SE = sqrt(Vgamma))
}





