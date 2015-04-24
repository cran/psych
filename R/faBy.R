"faBy" <- function(stats,nfactors=1,rotate="oblimin",fm="minres",free=TRUE,all=FALSE,min.n=12,...) {

 if(class(stats)[2] != "statsBy") stop("Please run statsBy first")
 cl <- match.call() 
  fo <- fa(stats$pooled,nfactors=nfactors,rotate=rotate,fm=fm) #get the overall pooled structure
  fl <- fo$loadings
  f <- list() #hold the results of each fa for each group
  ngroups <- stats$nG 
  nvar <- ncol(stats$r[[1]])
  
#replicateslist <- mclapply(1:ngroups,function(x) {
  
  stats$r <- pickgood(stats,min.n=min.n) 
  replicateslist <- lapply(stats$r,function(X,...) {
  if(!is.null(X) && !any(is.na(X))){
  fs <- try(fac(X,nfactors=nfactors,rotate=rotate,scores="none",...)) #call fa with the appropriate parameters
  if(class(fs)=="try-error") {warning("could not factor a within subject matrix")} else {
   if(!free && (nfactors > 1))  {fs$loadings <- faMatch(fl,fs)$loadings  #first match them
                                fs$loadings <- TargetQ(fs$loadings,Target=list(fl))$loadings} else {
                         fs$loadings <- faMatch(fl,fs)$loadings}
  if(!is.null(fs$Phi)) {  phis <- fs$Phi  
  if(all) { replicates <- list(fa=fs,loadings=(fs$loadings),phis=phis,vloadings = as.vector(fs$loadings),vphis = phis[lower.tri(phis)])}  else {
   replicates <- list(loadings=(fs$loadings),phis=phis,vloadings = as.vector(fs$loadings),vphis = phis[lower.tri(phis)])} } else 
     {replicates <- list(loadings=fs$loadings,vloadings <- as.vector(fs$loadings))}
  }} }) #end mclapply
  
  fabygroups <- lapply(replicateslist,function(X) X$vloadings)
  notnullgroup <- unlist(lapply(fabygroups,function(x) !is.null(x)))
  namesbygroup <- names(fabygroups)[notnullgroup] 
  fabygroups  <- matrix(unlist(lapply(replicateslist,function(X) X$vloadings)),ncol=nvar*nfactors,byrow=TRUE)
  num.groups <- nrow(fabygroups)
  means <- colMeans(fabygroups,na.rm=TRUE)
  sds <- apply(fabygroups,2,sd,na.rm=TRUE)
  faby <- matrix(c(as.vector(fl),means,sds,t(fabygroups)),ncol=(num.groups+3))
  fnames<- colnames(fo$loadings)[1:nfactors]
  vnames <- rownames(fo$loadings)

  colnames(faby) <- c("Total","Mean","sd",paste0("gr-",namesbygroup))
  rownames(faby) <- paste(rep(vnames,nfactors),"-",rep(fnames,each=nvar))
    if(!is.null(fo$Phi)) {
   vphis <- matrix(unlist(lapply(replicateslist,function(X) X$vphis)),nrow=num.groups,byrow=TRUE)
   means.phis <- colMeans(vphis)
   sds.phis <-     apply(vphis,2,sd,na.rm=TRUE) 
   phiby <- (matrix(c(fo$Phi[lower.tri(fo$Phi)],means.phis,sds.phis,t(vphis)), ncol=(num.groups+3),byrow=FALSE))
   
  
   colnames(phiby) <- c("Total","Mean","sd",paste0("gr-",namesbygroup))
   rownames(phiby) <-1:(nfactors*(nfactors-1)/2)
   k <- 1
   for (fi in 1:(nfactors-1)) { 
       for (fj in (fi+1):(nfactors)) {rownames(phiby)[k] <- paste(fnames[fi],"-",fnames[fj],sep="")
              k <- k +1 }}
              }
   if(all) {faBy <- list(fa=lapply(replicateslist,function(X) X$fa),loadings=faby,phis=phiby,Call=cl) } else {
    faBy <- list(loadings=t(faby),phis=t(phiby),Call=cl)}
    class(faBy) <- c("psych","faBy")
    return(faBy)
  }
  
  
  "faMatch" <- function(f1,f2) {
  fc <- factor.congruence(f1,f2)
  ord <- apply(fc,2,function(x) {which(x==max(x))})
  f2 <- fa.organize(f2,o=ord)
  return(f2)}
  
  "pickgood" <- function(stats,min.n) {
  new <- list()
 for (i in 1: length(stats$r)) {
    if(min(stats$n[[i]]) >= min.n) {new[i] <- stats$r[i]
    }   }
    names(new) <- names(stats$r) 
return(new)}