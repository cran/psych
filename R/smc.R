"smc" <-
function(R) {
 p <- dim(R)[2]
 if (dim(R)[1] != p) {R <- cor(R,use="pairwise")} else { R <- cov2cor(R)
      if (!is.matrix(R)) R <- as.matrix(R)}
 R.inv <- solve(R)
 smc <- 1 -1/diag(R.inv)
 return(smc)
 }