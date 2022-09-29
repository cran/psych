#probably not needed since the loadings of the fa.extend are functionally the beta weights
#redone to simply organize the results more clearly from fa.extend
#and to find the t and p values of the beta weights
#the beta weights are just factor extension  values

faReg <- faRegression <- function(r,nfactors=1,ov=NULL,dv=NULL, n.obs = NA, np.obs=NULL,correct=TRUE,rotate="oblimin",SMC=TRUE,warnings=TRUE, fm="minres",alpha=.1, omega=FALSE,cor="cor",use="pairwise",cor.correct=.5,weight=NULL,smooth=TRUE, ...) {
 cl <- match.call()
# if(is.numeric(ev)) ev  <- colnames(r)[ev]    #in case we are selecting variables 
 if(is.numeric(ov)) ov  <- colnames(r)[ov]
if(is.numeric(dv)) dv   <- colnames(r)[dv]
  nv <- c(ov, dv)
  nv  <- nfactors + length(dv)
  
  if(isCorrelation(r)) {
  #if the data are a correlation matrix, do an extension analysis
 #first find the extension analysis
# n.obs <- NROW(r) 
 fe <- fa.extend(r=r,nfactors=nfactors, ov = ov, ev = dv, n.obs = n.obs, np.obs= np.obs,
                correct=correct, rotate = rotate, SMC = SMC, warnings=warnings, fm =fm ,
                alpha = alpha, omega= FALSE,                                   
                cor="cor",use="pairwise",cor.correct=.5,weight=NULL,smooth=TRUE, ...) 
                 if(!isCorrelation(r))  {rdv <-  cor(r[,dv],use=use)} else { 
                                        rdv<- r[dv,dv]}
                                        fdv <- fe$loadings[dv,1:nfactors]
  } else {
  #do a conventional factor analysis on the ov and then find the factor score correlations with the dvs
  n.obs <- NROW(r)
  fe <- fa(r[ov], nfactors=nfactors, scores ="tenBerge")

  fdv <- cor(r[,dv],fe$scores,  use=use)
  rdv <- cor(r[,dv],use=use)
  }              
 Phi <- fe$Phi
#Phi <- diag(1,nrow=nfactors)  # we set this to identity to use the fa extension loadings

   #fdv <- fe$Structure[dv,1:nfactors]
 R <- matrix(NA,ncol=nv, nrow=nv)
 R[1:nfactors,1:nfactors] <- Phi
 R[1:nfactors,(nfactors +1):nv] <- t(fdv)
 R[(nfactors +1):nv, 1:nfactors]  <- (fdv)
 R[((nfactors +1):nv),((nfactors +1):nv)] <- rdv
rownames(R) <-  colnames(R) <- c(colnames(Phi) , dv)
 
 diag(R) <- 1
 set <- setCor(y=dv, x = colnames(fdv), data=R ,n.obs=n.obs, plot=FALSE)
 result <- list(regression=set,fa.extend =fe, dv.cor =rdv, R= R, Call=cl)
 class(result) <- c("psych", "fa.reg")
 return(result)
 #return(set) 
 }
 

  