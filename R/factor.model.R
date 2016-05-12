"factor.model" <-
function(f,Phi=NULL,U2=TRUE) {
   if(!is.matrix(f)) f <- as.matrix(f)
    if(is.null(Phi)) {Phi <- diag(1,dim(f)[2])}
    if(!is.matrix(Phi)) {Phi <- as.matrix(Phi)}
    if (!U2) diag(Phi) <- 1
    result<- f %*% Phi %*%  t(f)
     if (!U2) diag(result) <- 1
    return (result)}

