#modified Dec 10, 2008 to return 1 on diagonal if non-invertible
#modifed March 20, 2009 to return smcs * variance if covariance matrix is desired
#modified April 8, 2009 to remove bug introduced March 10 when using covar from data
#modified Jan 14, 2010 to test if matrix before cov2cor call.
#modified October 2, 2010 to convert smcs < 0 to 0 -- this is situation encountered with extreme missingness in sapa matrices
#modified April 23, 2015 to handle NAs in the correlation matrix
#smcs are found for the non-NA variables, then, smcs for the remaining ones are found from the  correlations for those with NAs
"smc" <-
function(R,covar =FALSE) {
failed=FALSE
wcc <- maxr <-  NULL
 p <- dim(R)[2]
 if(is.null(colnames(R))) colnames(R) <- paste0("V",1:p) 
 smc.all <- rep(NA,p)
 names(smc.all) <- colnames(R)
 if (dim(R)[1] != p) {if(covar) {C <- cov(R, use="pairwise")
                                 vari <- diag(C)
                                 R <- cov2cor(C)
                                 } else {R <- cor(R,use="pairwise")}}  else {vari <- diag(R)
                                 if (!is.matrix(R)) R <- as.matrix(R)
                                 R <- cov2cor(R)
                                 }
tempR <- NULL
if(any(is.na(R))) {
    bad <- TRUE
    tempR <- R
    vr <- diag(tempR)
    diag(tempR) <- 0

    maxr <- apply(tempR,2,function(x) max(abs(x),na.rm=TRUE))
    diag(tempR) <- vr 
    wcl <-NULL
    while(bad) {
     wc <- table(which(is.na(tempR), arr.ind=TRUE))  #find the correlations that are NA
    wcl <- c(wcl,as.numeric(names(which(wc==max(wc)))))
    tempR <- R[-wcl,-wcl]
    if(any(is.na(tempR))) {bad <- TRUE} else {bad <- FALSE}
         }
          warning("Missing values (NAs) in the correlation matrix do not allow for SMC's to be found for all variables.  \nI will try to estimate SMCs for those variables by their non-NA  correlations.")
     cat('\nSMC(s) for variables ',colnames(R)[wcl], 'were replaced (if possible) with smcs based upon their  (its) non-NA correlations\n')
     #now, try to find the smcs for the other ones
     wc <-(which(is.na(R[,wcl]),arr.ind=TRUE))

     if(is.null(dim(wc))) {wcc <- as.numeric(names(table(wc))) } else {
     wcc <- as.numeric(names(table(wc[,1])))}
     tempR <- R[-wcc,-wcc]  
      R <- R[-wcl,-wcl]
    
               }
if(!covar) { R <- cor.smooth(R) }  
                            
 R.inv <- try(solve(R),TRUE)
 if(class(R.inv)== as.character("try-error")) {smc <- rep(1,p)
 message("In smc, the correlation matrix was not invertible, smc's returned as 1s")} else  {smc <- 1 -1/diag(R.inv)}
 names(smc) <- colnames(R)
if(!is.null(tempR)) { R.na.inv <- try(solve(tempR),TRUE)
 if(class(R.na.inv) == as.character("try-error")) {smc.na <- rep(1,p)
   message("Oh bother, in smc, the correlation matrix of the adjusted part was not invertible, smc's returned as 1s")} else  {smc.na <- 1 -1/diag(R.na.inv)}
   } else {smc.na <- smc}
   
 if(all(is.na(smc))) {message ("Something is seriously wrong the correlation matrix.\nIn smc, smcs were set to 1.0")
                      smc[is.na(smc)] <- 1}
 if(max(smc,na.rm=TRUE) > 1.0) {message("In smc, smcs > 1 were set to 1.0")
   smc[smc >1 ]  <- 1.0}
   if(min(smc,na.rm=TRUE) < 0.0) {message("In smc, smcs < 0 were set to .0")
   smc[smc < 0]  <- 0}

   smc.all[names(smc.all) %in% names(smc)] <- smc

   if(!is.null(wcc)) {smc.all[wcl] <- smc.na[names(smc.all[wcl])] }
  
  smc <- smc.all
  
 if(!is.null(maxr)) { if(any(is.na(smc)))  {warning("The SMCs with NA values were replaced by their maximum correlation.") 
        cat('\nSMC(s) for variables ',names(smc)[is.na(smc)], 'were replaced with their maximum correlation \n')}
        smc[is.na(smc) ] <- maxr[is.na(smc)]  #in case we could not fix everything 
     }
 if(covar) {smc <- smc * vari}
 return(smc)
 }