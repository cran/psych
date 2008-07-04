"smc" <-
function(R) {
 p <- dim(R)[2]
 if (dim(R)[1] != p) {R <- cor(R,use="pairwise")} else { if (!is.matrix(R)) R <- as.matrix(R)}
 I <- diag(rep(1,p))
 R.inv <- solve(R)
 smc <- 1 -1/diag(R.inv)
 return(smc)
 }