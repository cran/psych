"irt.fa" <- 
function(x,nfactors=1,correct=TRUE,plot=TRUE,n.obs=NULL,...) {
cl <- match.call()
if (is.matrix(x) | is.data.frame(x)) {
	if(is.null(n.obs)) n.obs <- dim(x)[1]
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
t <- fa(r,nfactors=nfactors,n.obs=n.obs,...)
nf <- dim(t$loadings)[2]
 diffi <- list() 
     #flag <- which(abs(t$loadings) > 1,arr.ind=TRUE)
     #this throws an error if a Heywood case
     for (i in 1:nf) {diffi[[i]]  <- tau/sqrt(1-t$loadings[,i]^2)
     }
     
discrim <- t$loadings/sqrt(1-t$loadings^2)
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
"fa2irt" <- function(f,rho,plot=TRUE,n.obs=NULL) {
cl <- match.call()
tau <- rho$tau
r <- rho$rho
nf <- ncol(f$loadings)
diffi <- list() 
     #flag <- which(abs(t$loadings) > 1,arr.ind=TRUE)
     #this throws an error if a Heywood case
     for (i in 1:nf) {diffi[[i]]  <- tau/sqrt(1-f$loadings[,i]^2)
     }
discrim <- f$loadings/sqrt(1-f$loadings^2)
irt <- list(difficulty=diffi,discrimination=discrim)
result <- list(irt=irt,fa = f,rho=r,tau=tau,n.obs=n.obs,Call=cl)

if(class(rho)[2] == "poly" ) {class(result) <-  c("psych","irt.poly") } else {class(result) <- c("psych","irt.fa")}
if(plot) {pr <- plot(result) 
result$plot <- pr}
return(result)
}