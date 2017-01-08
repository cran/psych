"irt.fa" <- 
function(x,nfactors=1,correct=TRUE,plot=TRUE,n.obs=NULL,rotate="oblimin",fm="minres",sort=FALSE,...) {
cl <- match.call()
if (is.matrix(x) | is.data.frame(x)) {
	if(is.null(n.obs)) n.obs <- dim(x)[1]
	   nvar <- ncol(x)
	    vname <- colnames(x)
	 x <- as.matrix(x)
       if(!is.numeric(x)) {message("Converted non-numeric matrix input to numeric. \n Are you sure you wanted to do this?\n Please check your data")
         x <- matrix(as.numeric(x),ncol=nvar)}
             
	     colnames(x) <- vname 
	 
	tx <- table(as.matrix(x))
	if(dim(tx)[1] ==2) {tet <- tetrachoric(x,correct=correct)
	    typ = "tet"} else {tet <- polychoric(x)
	    typ = "poly"}

	r <- tet$rho
	tau <- tet$tau}  else {if (!is.null(x$rho)) { r <- x$rho
   			tau <- x$tau
   			if(is.null(n.obs)) {n.obs <- x$n.obs} 
             typ <- class(x)[2]
   			if (typ == "irt.fa") typ <- "tet"
   			 
   			  }  else {stop("x must  be a data.frame or matrix or the result from tetra or polychoric")}
              }
t <- fa(r,nfactors=nfactors,n.obs=n.obs,rotate=rotate,fm=fm,...)
if(sort) {t <- fa.sort(t)   #added 04/06/16
     if(typ !="tet" ) {tau <- tau[t$order,] } else {tau <- tau[t$order] }
     }

nf <- dim(t$loadings)[2]

 diffi <- list() 
     #flag <- which(abs(t$loadings) > 1,arr.ind=TRUE)
     #this throws an error if a Heywood case
     for (i in 1:nf) {diffi[[i]]  <- tau/sqrt(1-t$loadings[,i]^2)
     }
     
discrim <- t$loadings/sqrt(1-t$loadings^2)
if(any(is.nan(discrim))) {
      for (i in 1:nf) {
         bad <- which(is.nan(discrim[,i]))
        if(length( bad) > 0) { warning("A discrimination  with a NaN value was replaced with the maximum discrimination for factor ", i, " and  item(s)  ",bad, "\nexamine the factor analysis object (fa)  to identify the Heywood case. \nThe item informations are probably suspect as well for this factor.  \nYou might try a different factor extraction technique. ") 
          discrim[is.nan(discrim[,i]),i] <- max(discrim[,i],na.rm=TRUE) 
          diffi[[i]][bad,] <- tau[bad,]
         
         }
        }}
        

     


class(diffi) <- NULL
class(discrim) <- NULL
tl <- t$loadings
class(tl) <- NULL
irt <- list(difficulty=diffi,discrimination=discrim)
nlevels <- dim(diffi[[1]])[2]
#if(!is.null(nlevels)) {
#colnames(coeff) <- c(paste("Location",1:nlevels,sep=""),"Discrimination",paste("tau",1:nlevels,sep=""),"Loading") } else {
#colnames(coeff) <- c("Location","Discrimination","tau","Loading")}
result <- list(irt=irt,fa = t,rho=r,tau=tau,n.obs=n.obs,Call=cl)
switch(typ,
 tet = { class(result) <- c("psych","irt.fa")},
 tetra ={class(result) <- c("psych","irt.fa")},
 poly = {class(result) <- c("psych","irt.poly")},
 irt.poly = {class(result) <- c("psych","irt.poly")})

if(plot) {pr <- plot(result) 
result$plot <- pr}

return(result)
}

#convert a factor analysis output to an IRT output
#December 9, 2012
#modifed June 5, 2016 to allow sorted output to work
"fa2irt" <- function(f,rho,plot=TRUE,n.obs=NULL) {
cl <- match.call()
tau <- rho$tau
if(!is.null(f$order)) {if (!is.null(ncol(tau))) { tau <- tau[f$order,] } else {tau <- tau[f$order] }}  #adjust for sorting
r <- rho$rho
nf <- ncol(f$loadings)
diffi <- list() 
     #flag <- which(abs(t$loadings) > 1,arr.ind=TRUE)
     #this throws an error if a Heywood case
     for (i in 1:nf) {diffi[[i]]  <- tau/sqrt(1-f$loadings[,i]^2)
     }
discrim <- f$loadings/sqrt(1-f$loadings^2)
if(any(is.nan(discrim))) {bad <- which(is.nan(discrim),arr.ind=TRUE)
      if(length(bad) > 0) { warning("An discrimination  with a NaN value was replaced with the maximum discrimination for item(s)  ",bad, "\nexamine the factor analysis object (fa)  to identify the Heywood case") }
          for (i in 1:nf) {
          discrimin[is.nan(discrim[,i])] <- max(discrimin[,i],na.rm=TRUE) 
        }}
    
irt <- list(difficulty=diffi,discrimination=discrim)
result <- list(irt=irt,fa = f,rho=r,tau=tau,n.obs=n.obs,Call=cl)

if(class(rho)[2] == "poly" ) {class(result) <-  c("psych","irt.poly") } else {class(result) <- c("psych","irt.fa")}
if(plot) {pr <- plot(result) 
result$plot <- pr}
return(result)
}