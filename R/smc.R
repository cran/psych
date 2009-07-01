#modified Dec 10, 2008 to return 1 on diagonal if non-invertible
#modifed March 20, 2009 to return smcs * variance if covariance matrix is desired
#modified April 8, 2009 to remove bug introduced March 10 when using covar from data
"smc" <-
function(R,covar =FALSE) {
failed=FALSE
 p <- dim(R)[2]
 if (dim(R)[1] != p) {if(covar) {C <- cov(R, use="pairwise")
                                 vari <- diag(C)
                                 R <- cov2cor(C)
                                 } else {R <- cor(R,use="pairwise")}}  else {vari <- diag(R)
                                 R <- cov2cor(R)
      if (!is.matrix(R)) R <- as.matrix(R)}
 R.inv <- try(solve(R),TRUE)
 if(class(R.inv)== as.character("try-error")) {smc <- rep(1,p)
 message("In smc, the correlation matrix was not invertible, smc's returned as 1s")} else  {smc <- 1 -1/diag(R.inv)
 if(max(smc) > 1.0) {message("In smc, smcs > 1 were set to 1.0")
   smc[smc>1]  <- 1.0}
 if(covar) {smc <- smc * vari}}
 return(smc)
 }