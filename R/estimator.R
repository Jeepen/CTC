library(survival)

est <- function(time,status,treatmentstart,endofstudy,treatmentlength,window,data){
    time <- eval(substitute(time),data)
    status <- eval(substitute(status),data)
    treatmentstart <- eval(substitute(treatmentstart),data)
    expCases <- (status==1&(treatmentlength+window)<T)
    caseT <- time[expCases]
    caseL <- treatmentstart[expCases]

    casesExpRef <- ((caseT-window)>caseL&(caseT-window)<(caseL+treatmentlength))
    casesExpEvent <- (caseT>caseL&(caseT-treatmentlength)<caseL)
    nCases <- length(caseT)                           
    id <- rep(1:nCases,2)                              
    ## Collect in data.frame
    dCase <- data.frame(id=id,status=rep(c(1,0),each=nCases),
                        Exp=c(casesExpEvent,casesExpRef))
    tmp <- clogit(status~Exp+strata(id),data=dCase)$coefficients
    LCtrl <- treatmentstart[time>endofstudy]
    tmp - log((mean(LCtrl<endofstudy)-mean(LCtrl<treatmentlength))/
              mean(LCtrl<(endofstudy-treatmentlength)))
}
