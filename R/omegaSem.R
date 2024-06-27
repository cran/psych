#August 12, 2020  dropped the calls to the sem package.  Just use lavaan for the processing

"omegaSem" <-
function(m,nfactors=3,fm="minres",key=NULL,flip=TRUE, digits=2,title="Omega",
		sl=TRUE,labels=NULL, plot=TRUE,n.obs=NA,rotate="oblimin",
		Phi = NULL,option="equal",lavaan=TRUE,...) {
      #m is a correlation matrix, or if not, the correlation matrix is found
      #nfactors is the number of factors to extract
      #key allows items to be reversed scored  if desired
      #if Phi is not null, this implies that we have been given a factor matrix  -- added May 30, 2010
if(lavaan) {if(!requireNamespace('lavaan')) stop("You must have the lavaan package installed to use omegaSem")} #else {if(!requireNamespace('sem')) stop("You must have the sem package installed to use omegaSem")}
if (!sl) {warning("OmegaSem only works for Bifactor models, sl set to TRUE ")
    sl <- TRUE}

    cl <- match.call() 
om <- omega(m=m,nfactors=nfactors,fm=fm,key=key,flip=flip, digits=digits,
	title=title,sl=sl,labels=labels, plot=plot,n.obs=n.obs,
	rotate=rotate,Phi=Phi,option=option,...)
      #m is a correlation matrix, or if not, the correlation matrix is found
      #nfactors is the number of factors to extract
      #key allows items to be reversed scored  if desired
      #if Phi is not null, this implies that we have been given a factor matrix  -- added May 30, 2010
      
if(lavaan) {sem.model <- om$model$lavaan} else {sem.model <- om$model$sem}

if(is.na(n.obs)) {n.obs <- om$gstats$n.obs}

 if(dim(m)[1] != dim(m)[2]) {
                            n.obs <- dim(m)[1]
                            m <- cor(m,use="pairwise")} else {
                            m <- cov2cor(as.matrix(m))    #make sure it is a correlation matrix not a covariance or data matrix (if we change this, we will need to change the calculation for omega later)
                           }
      nvar <- dim(m)[2]
if(is.na(n.obs)) {message("Number of observations not specified. Arbitrarily set to 500")
                 n.obs <- 500
                 }
                            
      if(is.null(colnames(m))) {  rownames(m) <- colnames(m) <- paste("V",1:nvar,sep="") }
       m.names <- colnames(m)
if(lavaan) {if(!requireNamespace('lavaan')) stop("You must have the lavaan package installed to use omegaSem")} #else {if(!requireNamespace('sem')) stop("You must have the sem package installed to use omegaSem")}

#if(!requireNamespace('sem')) {stop("You must have the sem package installed to use omegaSem")
  
 if(lavaan) {sem.om <- lavaan::cfa(sem.model,sample.cov=m,sample.nobs=n.obs,orthogonal=TRUE,std.lv=TRUE) } #else {sem.om <- sem::sem(sem.model,m, n.obs) }
omega.efa <- omegaFromSem(sem.om,m,flip=flip,plot=plot)
results <- list(omegaSem=om,omega.efa=omega.efa,sem=sem.om,Call=cl)
class(results) <- c("psych", "omegaSem")
return(results)
}
#modified Sept 1 to pass the number of observations to SEM
#modified Jan 2, 2015 to call sem::sem  which seems to be the preferred manner

"omegaFromSem" <- 
function(fit,m=NULL,flip=TRUE,plot=TRUE) {
# m is the correlation matrix
# s is the sem solution from either sem or from lavaan -- I think I want to drop the sem option
if(inherits(fit,"lavaan")) { #get the lavaan parameters
  fx <- fit@Model@GLIST$lambda 
  m <- cov2cor(as.matrix(fit@SampleStats@cov[[1]]))
  Fit <- fit@Fit@test
  sem <- "lavaan"
  nvar <- nrow(fx)
  n.fact <- ncol(fx)
  g <- fx[,1]
  rn <- fit@Data@ov.names[[1]]
  cfa.loads <- fx} else {
  #get the sem parameters
  sem <- "sem"
nvar <- dim(m)[1]
n.fact <-dim(fit$A)[2]-nvar
rn <- rownames(m)
g <- fit$A[1:nvar,nvar+n.fact]
Fit <- NULL  #for now

cfa.loads <- cbind(g,fit$A[1:nvar,(nvar+1):(nvar+n.fact-1)])  #this puts g first
}

if(flip) { 
	flipper <- rep(1,nvar)
	flipper[g < 0] <- -1
	signkey <- strtrim(flipper,1)
            		 signkey[signkey=="1"] <- ""
            		 rn <- paste(rn,signkey,sep="")
	flipper <- diag(flipper)
	m <- flipper %*% m %*% t(flipper) 
	g <- flipper %*% g
	cfa.loads <- flipper %*% cfa.loads
	countneg <- colSums(cfa.loads < 0)
	for(i in 1:n.fact) {
	  if(countneg[i]> 0 ) cfa.loads[,i] <- -cfa.loads[,i] }
  }


w <- solve(m,cfa.loads)
rownames(cfa.loads)  <- rn
rownames(w)  <- rn
if(NCOL(cfa.loads) >1 ) {colnames(cfa.loads) <- c("gs",paste0("F",1:(n.fact-1),"s*"))} else {
warning("\nJust one factor defined, therefore, just omega_total is found")
 colnames(cfa.loads) <- "gs"}
gR2 <- diag(t(w) %*% cfa.loads)
Vt <- sum(m)
omh.sem <- sum(g)^2/Vt
h2 <- sum(rowSums(cfa.loads^2))
uniq <- tr(m) - h2
omt.sem <- (Vt - uniq)/Vt

 #find the subset omegas
     omg <- omgo <- omt<-  rep(NA,n.fact)
     sub <- apply(cfa.loads,1,function(x) which.max(abs(x[2:(n.fact)])))
     grs <- 0
 for(group in( 1:n.fact)) {
     groupi <- which(sub==group)
     if(length(groupi) > 0) {
      Vgr <- sum(m[groupi,groupi])
      gr <- sum(cfa.loads[groupi,(group+1)])
      grs <- grs + gr^2
      omg[group+1] <- gr^2/Vgr
      omgo[group+1] <- sum(cfa.loads[groupi,1])^2/Vgr
      omt[group+1] <- (gr^2+ sum(cfa.loads[groupi,1])^2)/Vgr
     }
     }
     omgo[1] <- sum(cfa.loads[,1])^2/sum(m)  #omega h
     omg[1] <- grs/sum(m)  #omega of subscales
     omt[1] <- omt.sem 
     om.group <- data.frame(total=omt,general=omgo,group=omg)


   #  rownames(om.group) <- colnames(gf$sl)[1:(nfactors+1)]
class(cfa.loads) <- "loadings"
results <- list(omega=omh.sem,omega.tot = omt.sem,cfa.loads=cfa.loads,gR2=gR2,omega.group=om.group,Fit=Fit,sem=sem)
class(results) <- c("psych","omegaSem")
if(plot) {if(n.fact > 1) { omega.diagram(results,sort=TRUE)} else {if(sem == 'lavaan') lavaan.diagram(fit,cut=0) } }
return(results)
}
#added the ability to work with just 1 factor (following error report by Erich Studerus) August 12, 2020



