#added checks for GPArotation even though we are already testing somewhere else
#the first function finds the first derivative, the second finds the fit
"vgbQ.bimin" <- function(L) {
k <- dim(L)[2]
L2 <- L^2
N <- matrix(1,k,k) 
diag(N) <- 0
L2N <- L2 %*% N
v <- sum(L2 * L2N)
G = 4 * L * L2N
return(list(f=v,Gq=G))
}

"vgQ.bimin" <- function(L) {
L2 <- L[,-1]
lvg <- vgbQ.bimin(L2)
v <- lvg$f
G <- lvg$Gq
G <-cbind(G[,1],G)
G[,1] <- 0
return(list(f=v,Gq=G))
}

#adapted from Jennrich and Bentler 2011
#requires GPArotation
"bifactor" <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000){
  if(requireNamespace('GPArotation')) {GPArotation::GPForth(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
           method="bimin")} else {stop("Bifactor requires GPArotation")}
   }
 #the oblique case  
#requires GPArotation
"biquartimin" <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000){
   if(requireNamespace('GPArotation')) {GPArotation::GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
           method="bimin") } else {stop("biquartimin requires GPArotation")}
   }
   

#this is a minor patch to the target function to allow it to have missing elements in the target so it more closely approximates the Michael Brown function  
"vgQ.targetQ" <- function (L, Target = NULL) 
{
    if (is.null(Target)) 
        stop("argument Target must be specified.")
       Gq <-  2 * (L - Target)
       Gq[is.na(Gq)] <- 0
    list(Gq = Gq, f = sum((L - Target)^2,na.rm=TRUE), Method = "Target rotation")
}

"TargetQ" <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000,Target=NULL) {
    if(requireNamespace('GPArotation')) {GPArotation::GPFoblq(L, Tmat=Tmat,normalize=normalize, eps=eps, maxit=maxit,
           method="targetQ",Target)} else {stop("TargetQ requires GPArotation")}}
           
"TargetT" <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000,Target=NULL) {
    if(requireNamespace('GPArotation')) {GPArotation::GPForth(L, Tmat=Tmat,normalize=normalize, eps=eps, maxit=maxit,
           method="targetQ",Target)} else {stop("TargetT requires GPArotation")}}

           
