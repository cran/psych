#modified Dec 10, 2008 to return 1 on diagonal if non-invertible
"smc" <-
function(R) {
failed=FALSE
 p <- dim(R)[2]
 if (dim(R)[1] != p) {R <- cor(R,use="pairwise")} else { R <- cov2cor(R)
      if (!is.matrix(R)) R <- as.matrix(R)}
 R.inv <- try(solve(R),TRUE)
 if(class(R.inv)== as.character("try-error")) {smc <- rep(1,p)
 warning("Correlation matrix not invertible, smc's returned as 1s")} else  {smc <- 1 -1/diag(R.inv)}
 return(smc)
 }