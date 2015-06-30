"faBy" <- function(stats,nfactors=1,rotate="oblimin",fm="minres",free=TRUE,all=FALSE,min.n=12,quant=.1,...) {

 if(class(stats)[2] != "statsBy") stop("Please run statsBy first")
 cl <- match.call() 
  fo.orth <- fa(stats$pooled,nfactors=nfactors,rotate="none",fm=fm) #get the overall pooled structure
  fo.rotated <-  fa(stats$pooled,nfactors=nfactors,rotate=rotate,fm=fm,...)  
  #could replace with a call to 
  #fo.rotated <- faRotate(fo.orth,rotate=rotate,...)
  # fl <- fo.rotated$loadings
    fl <- fo.rotated$loadings
  f <- list() #hold the results of each fa for each group
  ngroups <- stats$nG 
  nvar <- ncol(stats$r[[1]])
  
#replicateslist <- mclapply(1:ngroups,function(x) {
  
  
  stats$r <- pickgood(stats,min.n=min.n) #get the good analyses
  replicateslist <- lapply(stats$r,function(X,...) {
  

   if(!is.null(X) ){
  	 if(!free && (nfactors > 1)) {
         fs <- try(fac(X,nfactors=nfactors,rotate="none",scores="none",...))  #call fa but do not rotate
     	 #First match the orthogonal factors to get the right order and directions
       	#then target rotate the subject solution to the oblique pooled solution
      	#then match the order and directions (fixing the correlations) of this new solution
                                 
        fs$loadings <- faMatch(fo.orth,fs)$loadings  #first match them and then target rotate to overall
        fs<- TargetQ(fs$loadings,Target=list(fl))
        fs <- faMatch(fl,fs)  #try to match them but don't force a rotation
        }    else {
        fs <- try(fac(X,nfactors=nfactors,rotate=rotate,scores="none",...)) #call fa with the appropriate parameters
        fs$loadings <- faMatch(fl,fs)$loadings  #try to match them but don't force a rotation
      } 
  
  #if( length(class(fs)) ==1  ) {warning("could not factor a within subject matrix")} else {
  #     if(!free && (nfactors > 1))  { else {
   
    
  if(!is.null(fs$Phi)) {  phis <- fs$Phi 
 
  if(all) { replicates <- list(fa=fs,loadings=(fs$loadings),phis=phis,vloadings = as.vector(fs$loadings),vphis = phis[lower.tri(phis)])}  else {
  
   replicates <- list(loadings=fs$loadings,phis=phis,vloadings = as.vector(fs$loadings),vphis = phis[lower.tri(phis)])} 
   #} else 
   #  {replicates <- list(loadings=fs$loadings,vloadings <- as.vector(fs$loadings))}
 # }} 
 }
 }

 }
  ) #end mclapply
  

  fabygroups <- lapply(replicateslist,function(X) X$vloadings)
  notnullgroup <- unlist(lapply(fabygroups,function(x) !is.null(x)))
  namesbygroup <- names(fabygroups)[notnullgroup] 
  fabygroups  <- matrix(unlist(lapply(replicateslist,function(X) X$vloadings)),ncol=nvar*nfactors,byrow=TRUE)
  num.groups <- nrow(fabygroups)
  means <- colMeans(fabygroups,na.rm=TRUE)
  sds <- apply(fabygroups,2,sd,na.rm=TRUE)
  quants.low <- apply(fabygroups,2,quantile,quant)
  quants.high<- apply(fabygroups,2,quantile,1-quant)
  
  fnames<- colnames(fo.rotated$loadings)[1:nfactors]
  vnames <- rownames(fo.rotated$loadings)
 
  faby.sum <- matrix(c(as.vector(fl),means,sds,quants.low,quants.high),ncol=5)
  colnames(faby.sum) <-c("Pooled","mean","sd","low","high")
  rownames(faby.sum) <- paste(rep(vnames,nfactors))
  faby <- t(fabygroups)


  colnames(faby) <- c(paste0("gr-",namesbygroup))
  rownames(faby) <- paste(rep(vnames,nfactors),"-",rep(fnames,each=nvar))
    if(!is.null(fo.rotated$Phi)) {
   vphis <- matrix(unlist(lapply(replicateslist,function(X) X$vphis)),nrow=num.groups,byrow=TRUE)
   means.phis <- colMeans(vphis)
   sds.phis <-     apply(vphis,2,sd,na.rm=TRUE) 
   phis.low <- apply(vphis,2,quantile,quant)
   phis.high <- apply(vphis,2,quantile,1-quant)
   phiby.sum <- matrix(c(fo.rotated$Phi[lower.tri(fo.rotated$Phi)],means.phis,sds.phis,phis.low,phis.high),ncol=5)
   phiby <- (matrix(c(fo.rotated$Phi[lower.tri(fo.rotated$Phi)],means.phis,sds.phis,phis.low,phis.high,t(vphis)), ncol=(num.groups+5),byrow=FALSE))
   
  
   colnames(phiby) <- c("Total","Mean","sd","low","high", paste0("gr-",namesbygroup))
   rownames(phiby) <-1:(nfactors*(nfactors-1)/2)
   k <- 1
   for (fi in 1:(nfactors-1)) { 
       for (fj in (fi+1):(nfactors)) {rownames(phiby)[k] <- paste(fnames[fi],"-",fnames[fj],sep="")
              k <- k +1 }}
              }
    meanloading <- matrix(means,ncol=nfactors) 
    colnames(meanloading) <- fnames
    rownames(meanloading) <- vnames
    
     phis <- matrix(0,nfactors,nfactors)
     phis[lower.tri(phis)]   <- means.phis
     phis <-phis + t(phis)
     diag(phis) <- 1   
     colnames(phis) <- rownames(phis) <- fnames
   if(all) {faBy <- list(fa=lapply(replicateslist,function(X) X$fa),loadings=faby,Phi=phiby,Call=cl) } else {
    faBy <- list(mean.loading= meanloading,mean.Phi= phis,faby.sum=faby.sum,Phi.sum = phiby.sum,loadings=t(faby),Phi=t(phiby),nfactors=nfactors,quant,Call=cl)}
    class(faBy) <- c("psych","faBy")
    return(faBy)
  }
  
  
  "faMatch" <- function(f1,f2) {
  fc <- factor.congruence(f1,f2)
  ord <- 1:ncol(fc)
  for(i in 1:(ncol(fc)-1)) {
  new <- which(abs(fc[i,])==max(abs(fc[i,])) )
  old <- ord[i]
  ord[i] <- new
  ord[new] <- old }
  flip <- rep(1,ncol(fc))
  for (i in 1:ncol(fc)) {
  
  if(fc[ord[i],ord[i]] < 0) {                     
  f2$loadings[,ord[i]] <- f2$loadings[,ord[i]] * -1
  flip[i] <- -1
}
 if(!is.null(f2$Phi)) f2$Phi <- diag(flip) %*% f2$Phi %*% diag(flip)
 }
 # ord <- apply(fc,2,function(x) {which(abs(x)==max(abs(x)))})
  f2 <- fa.organize(f2,o=ord)
  return(f2)}
  
  "pickgood" <- function(stats,min.n) { #just look at those cases with good data
  new <- list()
 for (i in 1: length(stats$r)) {
    if(!any(is.na(stats$r[[i]])) & (min(stats$n[[i]]) >= min.n))  {new[i] <- stats$r[i]
    }  
   # if(min(stats$n[[i]]) >= min.n) 
     }
    names(new) <- names(stats$r) 
return(new)}