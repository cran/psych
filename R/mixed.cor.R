

#mixedCor was added April 28, 2017 to make mixed.cor easier to use
"mixed.cor" <- 
function(x=NULL,p=NULL,d=NULL,smooth=TRUE,correct=.5,global=TRUE,ncat=8,use="pairwise",method="pearson",weight=NULL)  {
cat("\nmixed.cor is deprecated, please use mixedCor.")
mixedCor(data=x,c=NULL,p=p,d=d,smooth=smooth,correct=correct,global=global,ncat=ncat,use=use,method=method,weight=weight)
}

"mixedCor" <- function(data=NULL,c=NULL,p=NULL,d=NULL,smooth=TRUE,correct=.5,global=TRUE,ncat=8,use="pairwise",method="pearson",weight=NULL) {
cl <- match.call() 
original <- colnames(data)
organize <- FALSE   #the default is to specify the continuous, the polytomous and the dichotomous
if((missing(c) | is.null(c)) && (missing(p) | is.null(p))  && (missing(d) | is.null(d))) { #figure out which kinds of variables we are using 
organize <- TRUE
nvar <- ncol(data)
data <- as.matrix(data)  #to speed things up
#first, check for data coded as factors,  these mess up.   We should quit with warning 
#also check for variables with no variance and flag them (don't remove, just quit.  Likely problem with data)

progressBar(nvar,nvar,"Preparing the data")
ans <- matrix(NA,nrow=nvar,ncol=2)
   for (i in 1:nvar) {
        if (is.numeric(data[,i])) {ans[i,2] <- 1 } else {
            if ((is.factor(data[,i])) || (is.logical(data[,i]))) {
               ans[i,2]  <- 2
            } else {
                if (is.character(data[,i])) {
               ans[i,2] <- 3
                } else {ans[i,2] <- 4}
            }
        }
       ans[i,1] <- sd(data[,i],na.rm=TRUE)
    }
    if(any(ans[,2] !=1)) {cat("\nNot all of the variables are numeric.  Please check your data. \nPotential bad items are ")
    print(which(ans[,2] !=1)) 
    stop ("\nI am stopping because of the problem with the data.")
    }
    
       
 bad <- which(ans[,1] ==0)   
 if(length(bad) > 0) {cat("\nSome of the variables have no variance.  Please remove them and try again.\nBad items are ")
      print(bad)
      stop("\nI am stopping because of the problems with the data.")
      }


tab <- apply(data,2,table)
if(is.list(tab)) {len <- lapply(tab,length)} else {len <- dim(tab)[1] }
dvars <- subset(1:nvar,len==2)   #find the dichotomous variable by number
pvars <- subset(1:nvar,((len > 2) & (len <= ncat)))  #find the polytomous variables
cvars <- subset(1:nvar,(len > ncat))  #find the continuous variables (more than ncat levels)

#This next part is not as efficient as it should be, because it is making up 3 new matrices, where we could just do it by reference -- pass the names, not the data
#if(length(dvars) > 0) {d <- as.matrix(x[,dvars],ncol=length(dvars))
#              colnames(d) <- colnames(x)[dvars]} else {d <- NULL}
if(length(pvars) > 0) {#p <- as.matrix(x[,pvars],ncol=length(pvars))
              # colnames(p) <- colnames(x)[pvars] 
              # tab <- table(p) #now check to make sure that they are all on the same scale
              # if(length(tab) > ncat) stop("I tried to figure out which were continuous and which were polytomous, but failed.  Please try again by specifying x, p, and d.")
              tab <- apply(data[,pvars,drop=FALSE],2,table)  #if all elements have the same number of categories, this is a matrix, otherwise, it is a list
             if(is.list(tab)) {
               ok <- apply(data[,pvars,drop=FALSE], 2,function (x) {if (length(table(x)) != (max(x,na.rm=TRUE) - min(x,na.rm=TRUE)+1)) {FALSE} else {TRUE}})
               if(any(!ok)) {bad <- which(!ok)
               cat("\n Some polytomous variables have fewer categories than they should.  Please check your data.  \nPotential bad items are ",colnames(p)[bad],"\n")
             
               stop("\nI am stopping because of the problem with polytomous data")
               }
               }
                         } else {p <- NULL}
# if(length(cvars) > 0) {cont <- matrix(x[,cvars],ncol=length(cvars))
#                       colnames(cont) <- colnames(x)[cvars] } else {cont <- NULL}
                       
#Rho <- mixed.cor1(cont,p, d,smooth=smooth,global=global,correct=correct,use=use,method=method,weight=weight)
progressBar(nvar,nvar,"Starting mixed.cor1")
Rho <- mixed.cor1(data,cvars,pvars, dvars,smooth=smooth,global=global,correct=correct,use=use,method=method,weight=weight)
oldorder <- c(cvars,pvars,dvars)
ord <- order(oldorder)
Rho$rho <- Rho$rho[ord,ord]
} else {# organization is specified
#if ((p+ d) == nvar)
#if(!missing(c)) x <- data[c]
#if(!missing(p)) p <- data[p]
#if(!missing(d)) d <- data[d]
progressBar(1,1,"Starting mixed.cor1")
Rho <- mixed.cor1(data,c=c,p=p,d=d,smooth=smooth,global=global,correct=correct,use=use,method=method,weight=weight)}
orig <- original %in%  colnames(Rho$rho)
orig <- original[orig]
Rho$rho <- Rho$rho[orig,orig]  #organize them in the way they came in
Rho$Call <- cl
return(Rho)
}

#modified April 29th (2014) to get around the problem of missing rows or columns in the polychoric function.
#modified 1/1/14 to add multicore capability
#deprecated as of 02/05/18
# "mixed.cor" <- 
# function(x=NULL,p=NULL,d=NULL,smooth=TRUE,correct=.5,global=TRUE,ncat=8,use="pairwise",method="pearson",weight=NULL)  {
# cl <- match.call() 
# organize <- FALSE   #the default is to specify the continuous, the polytomous and the dichotomous
# if(!is.null(x) && is.null(p)  && is.null(d)) { #figure out which kinds of variables we are using 
# organize <- TRUE
# nvar <- ncol(x)
# x <- as.matrix(x)
# tab <- apply(x,2,function(x) table(x))
# if(is.list(tab)) {len <- lapply(tab,function(x) length(x))} else {len <- dim(tab)[1] }
# dvars <- subset(1:nvar,len==2)   #find the dichotomous variables
# pvars <- subset(1:nvar,((len > 2) & (len <= ncat)))  #find the polytomous variables
# cvars <- subset(1:nvar,(len > ncat))  #find the continuous variables (more than ncat levels)
# 
# if(length(dvars) > 0) {d <- as.matrix(x[,dvars],ncol=length(dvars))
#               colnames(d) <- colnames(x)[dvars]} else {d <- NULL}
# if(length(pvars) > 0) {p <- as.matrix(x[,pvars],ncol=length(pvars))
#                colnames(p) <- colnames(x)[pvars] 
#                tab <- table(p) #now check to make sure that they are all on the same scale
#                if(length(tab) > ncat) stop("I tried to figure out which were continuous and which were polytomous, but failed.  Please try again by specifying x, p, and d.")
#                ok <- apply(p, 2,function (x) {if (length(table(x)) != (max(x,na.rm=TRUE) - min(x,na.rm=TRUE)+1)) {FALSE} else {TRUE}})
#                if(any(!ok)) {bad <- which(!ok)
#                cat("\n Some polytomous variables have fewer categories than they should.  Please check your data.  \nPotential bad items are ",colnames(p)[bad],"\n")
#              
#                stop("\nI am stopping because of the problem with polytomous data")
#                }
#                          } else {p <- NULL}
#  if(length(cvars) > 0) {cont <- matrix(x[,cvars],ncol=length(cvars))tab <- apply(x,2,table)
# if(is.list(tab)) {len <- lapply(tab,length)} else {len <- dim(tab)[1] }
# dvars <- subset(1:nvar,len==2)   #find the dichotomous variable by number
# pvars <- subset(1:nvar,((len > 2) & (len <= ncat)))  #find the polytomous variables
# cvars <- subset(1:nvar,(len > ncat))  #find the continuous variables (more than ncat levels)

#                        colnames(cont) <- colnames(x)[cvars] } else {cont <- NULL}
#                        
# Rho <- mixed.cor1(cont,p, d,smooth=smooth,global=global,correct=correct,use=use,method=method,weight=weight)
# oldorder <- c(cvars,pvars,dvars)
# ord <- order(oldorder)
# Rho$rho <- Rho$rho[ord,ord]
# } else {# organization is specified
# #if ((p+ d) == nvar) 
# Rho <- mixed.cor1(x=x,p=p,d=d,smooth=smooth,global=global,correct=correct,use=use,method=method,weight=weight)}
# 
# Rho$Call <- cl
# return(Rho)
# }

#December 22,2010
#revised July 15, 2011 to work for the various special cases
#meant to combine continuous, polytomous and dichotomous correlations
#revised October 12, 2011 to get around the sd of vectors problem
#revised Sept 10, 2013 to allow for dichotomies with different minima	
"mixed.cor1" <-
function(data,c=NULL,p=NULL,d=NULL,smooth=TRUE,global=TRUE,correct=correct,use=use,method=method,weight=NULL) {
 cl <- match.call() 
 x <- c   #the continuous variables
# if(!is.null(x)) {nx <- dim(x)[2]} else {nx <- 0}
# if(!is.null(p)) {np <- dim(p)[2]} else {np <- 0}
# if(!is.null(d))  {nd <- dim(d)[2]} else {nd <- 0}
# if(is.null(nx)) nx <- 1
# if(is.null(np)) np <- 1
# if(is.null(nd)) nd <- 1
if(!is.null(x)) {nx <- length(x)} else {nx <- 0}
if(!is.null(p)) {np <- length(p)} else {np <- 0}
 if(!is.null(d))  {nd <-length(d)} else {nd <- 0}
npd  <- nx +np + nd
rpd <- NULL #in case we don't have any
#check to make sure all the data are ok for doing the appropriate analyses
#first check for polychorics
#but we did this before in mixedCor.  If the person specified  the p and d, assume they are correct
# if(np > 0) { ptable <- table(as.matrix(data[,p]))   #haven't we already done this?
# nvalues <- length(ptable)  #find the number of response alternatives 
# if(nvalues > 8) stop("You have more than 8 categories for your items, polychoric is probably not needed")}
# 
# #now test for tetrachorics
# if(nd > 0) {
# dm <- apply(data[,d,drop=FALSE],2,function(x) min(x,na.rm=TRUE))
# data[,d]  <- t(t(data[,d]) - dm)  
# #d <- d -min(d,na.rm=TRUE) #in case the numbers are not 0,1   This requires them all to be the same
# if(max(data[,d,drop=FALSE],na.rm=TRUE) > 1) {stop("Tetrachoric correlations require dichotomous data")}}

if(nx > 1) {
progressBar(nx,nx,"Finding Pearson correlations")
rx <- cor(data[,x],use=use,method=method)} else {if(nx  < 1) {rx <- NULL
   rho <- NULL} else {rx <- 1
     rho <- 1}}
     
if(np > 1) {
#cat("\n Starting to find the polychoric correlations")
progressBar(np,np,"Finding polychorics")
 rp <- polychoric(data[,p],smooth=smooth,global=global,weight=weight,correct=correct)}    else {if (np == 1) {
	rho <- 1
	names(rho) <- colnames(p)
	rp <- list(rho=rho,tau=NULL)}  else {rp <- list(rho= NULL,tau=NULL)}}
if(nd > 1) {
# cat("n Starting to find the tetrachoric correlations\n\n")
progressBar(nd,nd,"Finding tetrachorics")
  rd <- tetrachoric(data[,d],smooth=smooth,correct=correct,weight=weight)}   else {if (nd == 1) {rd <- list(rho=1,tau=NULL)}  else {rd <- list(rho=NULL,tau=NULL)}}

if(nx > 0) {if(np > 0) {rxp <- polyserial(data[,x,drop=FALSE],data[,p,drop=FALSE])   #the normal case is for all three to exist
		tmixed <- cbind(rx,t(rxp))
		lmixed <- cbind(rxp,rp$rho)
		rho <- rbind(tmixed,lmixed)} else {rho <- rx}  #we now have x and p

		if(nd > 0) { rxd <- biserial(data[,x,drop=FALSE],data[,d,drop=FALSE])
		
		    if(np > 0) {
		    rpd <- polydi(data[,p,drop=FALSE],data[,d,drop=FALSE],global=global,correct=correct)$rho   #passing the global value (July 10, 2017)
			            topright <- t(cbind(rxd,t(rpd))) 
			            } else {
			        topright <- t(rxd)}
			tmixed <- cbind(rho,topright) 
			lmixed <- cbind(t(topright),rd$rho)
			rho <- rbind(tmixed,lmixed) }

		} else {  #the case of nx =0
		   if( np > 0) { 
		      if (nd > 0 ) {
		     progressBar(nd,nd,"Starting polydi")
		      rpd <- polydi(data[,p,drop=FALSE],data[,d,drop=FALSE],global=global,correct=correct)$rho  #added global  (July 10, 2017)   But this causes problems
		     
		       tmixed <- cbind(rp$rho,rpd)
		       lmixed <- cbind(t(rpd),rd$rho)
		       rho <- rbind(tmixed,lmixed)
		       }  else {rho <- rp$rho} } else {
		    rho <- rd$rho}
		    }
colnames(rho) <- rownames(rho)
class(rho) <- c("psych","mixed")
if(!is.null(rx)) class(rx) <- c("psych","mixed")
mixed <- list(rho=rho,rx=rx,poly=rp,tetra=rd,rpd=rpd,Call=cl)
class(mixed) <- c("psych","mixed")
return(mixed)
}

