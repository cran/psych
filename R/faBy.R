"faBy" <- function(stats,nfactors=1,rotate="oblimin",fm="minres",free=TRUE,all=FALSE,...) {
 if(class(stats)[2] != "statsBy") stop("Please run statsBy first")
 cl <- match.call() 
  fo <- fa(stats$pooled,nfactors=nfactors,rotate=rotate,fm=fm) #get the overall pooled structure
  fl <- fo$loadings
  f <- list() #hold the results of each fa for each group
  ngroups <- stats$nG 
  nvar <- ncol(stats$r[[1]])
  
#replicateslist <- mclapply(1:ngroups,function(x) {
  replicateslist <- lapply(stats$r,function(X,...) {
  fs <- fac(X,nfactors=nfactors,rotate=rotate,...) #call fa with the appropriate parameters
   if(!free && (nfactors > 1))  {fs$loadings <- target.rot(fs$loadings,fl)$loadings}
  if(!is.null(fs$Phi)) {  phis <- fs$Phi  
  if(all) { replicates <- list(fa=fs,loadings=(fs$loadings),phis=phis,vloadings = as.vector(fs$loadings),vphis = phis[lower.tri(phis)])}  else {
   replicates <- list(loadings=(fs$loadings),phis=phis,vloadings = as.vector(fs$loadings),vphis = phis[lower.tri(phis)])} } else 
     {replicates <- list(loadings=fs$loadings,vloadings <- as.vector(fs$loadings))}
  } ) #end mclapply
  
  fabygroups  <- matrix(unlist(lapply(replicateslist,function(X) X$vloadings)),nrow=ngroups,byrow=TRUE)
  means <- colMeans(fabygroups,na.rm=TRUE)
  sds <- apply(fabygroups,2,sd,na.rm=TRUE)
  faby <- matrix(c(as.vector(fl),means,sds,t(fabygroups)),ncol=(ngroups+3))
  fnames<- colnames(fo$loadings)[1:nfactors]
  vnames <- rownames(fo$loadings)
  colnames(faby) <- c("Total","Mean","sd",paste("gr-",1:ngroups,sep=""))
  rownames(faby) <- paste(rep(vnames,nfactors),"-",rep(fnames,each=nvar))
    if(!is.null(replicateslist[[1]]$vphis)) {
   vphis <- matrix(unlist(lapply(replicateslist,function(X) X$vphis)),nrow=ngroups,byrow=TRUE)
   means.phis <- colMeans(vphis)
   sds.phis <-     apply(vphis,2,sd,na.rm=TRUE) 
   phiby <- (matrix(c(fo$Phi[lower.tri(fo$Phi)],means.phis,sds.phis,t(vphis)), ncol=(ngroups+3),byrow=FALSE))
   colnames(phiby) <- c("Total","Mean","sd",paste("gr-",1:ngroups,sep=""))
   rownames(phiby) <-1:(nfactors*(nfactors-1)/2)
   k <- 1
   for (fi in 1:(nfactors-1)) { 
       for (fj in (fi+1):(nfactors)) {rownames(phiby)[k] <- paste(fnames[fi],"-",fnames[fj],sep="")
              k <- k +1 }}
              }
   if(all) {faBy <- list(fa=lapply(replicateslist,function(X) X$fa),loadings=faby,phis=phiby,Call=cl) } else {
    faBy <- list(loadings=faby,phis=phiby,Call=cl)}
    class(faBy) <- c("psych","faBy")
    return(faBy)
  }
  
  
  