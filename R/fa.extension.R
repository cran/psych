"fa.extension" <-
  function(Roe,fo,correct=TRUE) {
 cl <- match.call()
 omega <-FALSE

 if(!is.null(class(fo)[2])) {if(class(fo)[2]=="fa") {
          if(!is.null(fo$Phi)) {Phi <- fo$Phi} else {Phi <- NULL}
          
       fl <- fo$loadings 
       
     } else {if (class(fo)[2] == "omega") {
         omega <- TRUE
         w <- fo$stats$weights
         fl <- fo$schmid$sl
         Phi <- NULL
         fl <- fl[,1:(dim(fl)[2]-3)]
         nfactors <- dim(fl)[2]
          fe <- t(t(w) %*% Roe)
       } 
    }
    }
 
 if(!omega) fe <- t( Roe) %*% fl %*% (solve(t(fl)%*% (fl))) 
  if(!is.null(Phi)) fe <- fe %*% solve(Phi)
 
 if(!correct) {#the Gorsuch case
     d <-diag(t(fl) %*% fo$weight)
     fe <- (fe * d)
 }
 colnames(fe) <- colnames(fl)
if(!is.null(Phi)) {resid <- Roe - fl %*% Phi %*% t(fe)} else {resid <- fl  %*% t(fe)}
 result <- list(loadings = fe,Phi=Phi,resid=resid,Call=cl)
 class(result) <- c("psych","extension")
 return(result)
}
#written April 5, 2011
#revised August 15, 2011 to avoid using the weights matrix except in the omega case

