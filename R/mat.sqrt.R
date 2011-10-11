"mat.sqrt" <- function(x) {
   e <- eigen(x)
   sqrt.ev <- sqrt(e$values)   #need to put in a check here for postive semi definite
   inv.evec <- solve(e$vectors)
   result <- e$vectors %*% diag(sqrt.ev) %*% inv.evec
   result}

   