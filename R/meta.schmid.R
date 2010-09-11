"meta.schmid" <-
function(L = NULL,L2=NULL,g3 =NULL,m=NULL) {
Vt = sum(m)
g <- L %*% L2
omega <- sum(g)^2/Vt
result <- list(g,omega)
return(result)
}