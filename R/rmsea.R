RMSEA <- function(chisq, dof,n.obs,alpha=.1){
 tail <- alpha/2
if(!is.na(chisq)){
 RMSEA <- sqrt(max(chisq/(dof* n.obs) - 1/(n.obs-1), 0)) 

 max <- max(n.obs,chisq) +2* n.obs       
  RMSEA.U <- 0  #in case we can not find it
if(pchisq(df=dof,q=chisq) > tail){ RMSEA.U <-try ( sqrt(uniroot(function(x)    {pchisq(df=dof,ncp=x,q=chisq)- tail},c(0,max))$root/(n.obs-1)/dof),silent=TRUE)
    
if(inherits( RMSEA.U,"try-error")) {if(RMSEA <= 0 ) {RMSEA.U <- 0} else {message(" I could not find the RMSEA upper bound . Sorry about that")
       #if the fit is super good, then the chisq is too small to get an upper bound.  Report it as 0.
 RMSEA.U <- NA}}
         
    }

RMSEA.L <- 0  #in case we can not find it
   if(pchisq(df=dof,q=chisq) > (1-tail)) { RMSEA.L   <- try( sqrt(uniroot(function(x) {pchisq(df=dof,ncp=x,q=chisq)-1 + tail},c(0,max))$root/(n.obs-1)/dof) ,silent=TRUE)

if(inherits(RMSEA.L,"try-error")) { RMSEA.L <- NA} } else {RMSEA.L <- 0}
# 

result <- list(RMSEA.L=RMSEA.L,RMSEA=RMSEA,RMSEA.U=RMSEA.U,alpha=alpha) 


class(result) <- c("psych","rmsea")
return(result)    
} else {result <- RMSEA <- RMSEA.L <- RMSEA.U <- NULL}         
}

# chisq <- objective * ((n.obs-1) -(2 * n.var + 5)/6 -(2*nfactors)/3)