#modified Dec 10, 2008 to return 1 on diagonal if non-invertible
#modifed March 20, 2009 to return smcs * variance if covariance matrix is desired
#modified April 8, 2009 to remove bug introduced March 10 when using covar from data
#modified Jan 14, 2010 to test if matrix before cov2cor call.
#modified October 2, 2010 to convert smcs < 0 to 0 -- this is situation encountered with extreme missingness in sapa matrices
"smc" <-
function(R,covar =FALSE) {
failed=FALSE
 p <- dim(R)[2]
 if (dim(R)[1] != p) {if(covar) {C <- cov(R, use="pairwise")
                                 vari <- diag(C)
                                 R <- cov2cor(C)
                                 } else {R <- cor(R,use="pairwise")}}  else {vari <- diag(R)
                                 if (!is.matrix(R)) R <- as.matrix(R)
                                 R <- cov2cor(R)
                                 }
if(!covar) { R <- cor.smooth(R) }                              
 R.inv <- try(solve(R),TRUE)
 if(class(R.inv)== as.character("try-error")) {smc <- rep(1,p)
 message("In smc, the correlation matrix was not invertible, smc's returned as 1s")} else  {smc <- 1 -1/diag(R.inv)
 if(all(is.na(smc))) {message ("Something is seriously wrong the correlation matrix.\nIn smc, smcs were set to 1.0")
                      smc[is.na(smc)] <- 1}
 if(max(smc,na.rm=TRUE) > 1.0) {message("In smc, smcs > 1 were set to 1.0")
   smc[smc >1 ]  <- 1.0}
   if(min(smc,na.rm=TRUE) < 0.0) {message("In smc, smcs < 0 were set to .0")
   smc[smc < 0]  <- 0}
 if(covar) {smc <- smc * vari}}
 return(smc)
 }