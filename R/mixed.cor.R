#modified April 29th to get around the problem of missing rows or columns in the polychoric function.
#modified 1/1/14 to add multicore capability
"mixed.cor" <- 
function(x=NULL,p=NULL,d=NULL,smooth=TRUE,correct=.5,global=TRUE,ncat=8,use="pairwise",method="pearson",weight=NULL)  {
cl <- match.call() 
organize <- FALSE   #the default is to specify the continuous, the polytomous and the dichotomous
if(!is.null(x) && is.null(p)  && is.null(d)) { #figure out which kinds of variables we are using 
organize <- TRUE
nvar <- ncol(x)
x <- as.matrix(x)
tab <- apply(x,2,function(x) table(x))
if(is.list(tab)) {len <- lapply(tab,function(x) length(x))} else {len <- dim(tab)[1] }
dvars <- subset(1:nvar,len==2)   #find the dichotomous variables
pvars <- subset(1:nvar,((len > 2) & (len <= ncat)))  #find the polytomous variables
cvars <- subset(1:nvar,(len > ncat))  #find the continuous variables (more than ncat levels)

if(length(dvars) > 0) {d <- matrix(x[,dvars],ncol=length(dvars))
              colnames(d) <- colnames(x)[dvars]} else {d <- NULL}
if(length(pvars) > 0) {p <- matrix(x[,pvars],ncol=length(pvars))
               colnames(p) <- colnames(x)[pvars] 
               tab <- table(p) #now check to make sure that they are all on the same scale
               if(length(tab) > ncat) stop("I tried to figure out which where continuous and which were polytomous, but failed.  Please try again by specifying x, p, and d.")
               ok <- apply(p, 2,function (x) {if (length(table(x)) != (max(x,na.rm=TRUE) - min(x,na.rm=TRUE)+1)) {FALSE} else {TRUE}})
               if(any(!ok)) {bad <- which(!ok)
               cat("\n Some polytomous variables have fewer categories than they should.  Please check your data.  \nPotential bad items are ",colnames(p)[bad],"\n")
             
               stop("\nI am stopping because of the problem with polytomous data")
               }
                         } else {p <- NULL}
 if(length(cvars) > 0) {cont <- matrix(x[,cvars],ncol=length(cvars))
                       colnames(cont) <- colnames(x)[cvars] } else {cont <- NULL}
                       
Rho <- mixed.cor1(cont,p, d,smooth=smooth,global=global,correct=correct,use=use,method=method,weight=weight)
oldorder <- c(cvars,pvars,dvars)
ord <- order(oldorder)
Rho$rho <- Rho$rho[ord,ord]
} else {# organization is specified
#if ((p+ d) == nvar) 
Rho <- mixed.cor1(x=x,p=p,d=d,smooth=smooth,global=global,correct=correct,use=use,method=method,weight=weight)}

Rho$Call <- cl
return(Rho)
}

#December 22,2010
#revised July 15, 2011 to work for the various special cases
#meant to combine continuous, polytomous and dichotomous correlations
#revised October 12, 2011 to get around the sd of vectors problem
#revised Sept 10, 2013 to allow for dichotomies with different minima
"mixed.cor1" <-
function(x=NULL,p=NULL,d=NULL,smooth=TRUE,global=TRUE,correct=.5,use=use,method=method,weight=NULL) {
 cl <- match.call() 
if(!is.null(x)) {nx <- dim(x)[2]} else {nx <- 0}
if(!is.null(p)) {np <- dim(p)[2]} else {np <- 0}
if(!is.null(d))  {nd <- dim(d)[2]} else {nd <- 0}
if(is.null(nx)) nx <- 1
if(is.null(np)) np <- 1
if(is.null(nd)) nd <- 1
npd  <- nx +np + nd
#check to make sure all the data are ok for doing the appropriate analyses
#first check for polychorics
if(np > 0) { ptable <- table(as.matrix(p))
nvalues <- length(ptable)  #find the number of response alternatives 
if(nvalues > 8) stop("You have more than 8 categories for your items, polychoric is probably not needed")}

#now test for tetrachorics
if(nd > 0) {
dm <- apply(d,2,function(x) min(x,na.rm=TRUE))
d <- t(t(d) - dm)  
#d <- d -min(d,na.rm=TRUE) #in case the numbers are not 0,1   This requires them all to be the same
if(max(d,na.rm=TRUE) > 1) {stop("Tetrachoric correlations require dichotomous data")}}


if(nx > 0) {rx <- cor(x,use=use,method=method)} else {rx <- NULL
   rho <- NULL}
if(np > 1) {
#cat("\n Starting to find the polychoric correlations")
 rp <- polychoric(p,smooth=smooth,global=global,weight=weight)}    else {if (np == 1) {
	rho <- 1
	names(rho) <- colnames(p)
	rp <- list(rho=rho,tau=NULL)}  else {rp <- list(rho= NULL,tau=NULL)}}
if(nd > 1) {
 #cat("n Starting to find the tetrachoric correlations\n\n")
  rd <- tetrachoric(d,smooth=smooth,correct=correct,weight=weight)}   else {if (nd == 1) {rd <- list(rho=1,tau=NULL)}  else {rd <- list(rho=NULL,tau=NULL)}}

if(nx > 0) {if(np > 0) {rxp <- polyserial(x,p)   #the normal case is for all three to exist
		tmixed <- cbind(rx,t(rxp))
		lmixed <- cbind(rxp,rp$rho)
		rho <- rbind(tmixed,lmixed)} else {rho <- rx}  #we now have x and p
		
		if(nd > 0) { rxd <- biserial(x,d)
		    if(np > 0) {rpd <- polydi(p,d,correct=correct)$rho  
			            topright <- t(cbind(rxd,t(rpd))) 
			            } else {
			        topright <- t(rxd)}
			tmixed <- cbind(rho,topright) 
			lmixed <- cbind(t(topright),rd$rho)
			rho <- rbind(tmixed,lmixed) }

		} else {  #the case of nx =0
		   if( np > 0) { 
		      if (nd >0 ) {rpd <- polydi(p,d,correct=correct)$rho
		     
		       tmixed <- cbind(rp$rho,rpd)
		       lmixed <- cbind(t(rpd),rd$rho)
		       rho <- rbind(tmixed,lmixed)
		       }  else {rho <- rp$rho} } else {
		    rho <- rd$rho}
		    }
colnames(rho) <- rownames(rho)
class(rho) <- c("psych","mixed")
if(!is.null(rx)) class(rx) <- c("psych","mixed")
mixed <- list(rho=rho,rx=rx,poly=rp,tetra=rd,Call=cl)
class(mixed) <- c("psych","mixed")
return(mixed)
}

