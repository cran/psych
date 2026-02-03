#created 1/20/2026 as helper function for CFA  (code from omega -- might use in omega)
#need to handle a gmodel 
omegaStats <- function(r,loads =NULL,Phi=NULL,gload=NULL,group=NULL, h2=NULL,covar=FALSE,lavaan=NULL){
 	cl <- match.call()
 nvar <- ncol(r)
 if(!is.null(lavaan)) { fit <- lavaan     #grab the lavaan output
     loads <- fit@Model@GLIST$lambda
      colnames(loads) <-  fit@Model@dimNames[[3]][[1]]
      rownames(loads) <-  fit@Model@dimNames[[1]][[1]]
      Phi <- fit@Model@GLIST$psi
      rownames(Phi) <- colnames(Phi) <- colnames(loads)
      }
      ## normal case
 if(is.null(gload)) gload <- loads[,1]
 if(is.null(group)) group <- loads[,-1] 
 if(is.null(loads) ) loads <- cbind(gload,group)
 if(is.null(h2))  h2 <- rowSums(loads^2)   
# if(is.null(Phi)) {h2 <- diag(loads %*% t(loads) )} else {h2 <- diag(group %*% Phi %*% t(group))}

 nfactors <- ncol(group)
 Vt <- sum(r)   #find the total variance in the scale
 Vitem <- sum(diag(r)) 
  gsq <- (sum(gload))^2
   u2 <- diag(r) - h2
    uniq <- sum(u2)
     om.tot <- (Vt-uniq)/Vt
      om.limit <- gsq/(Vt-uniq)
     alpha <- ((Vt-Vitem)/Vt)*(nvar/(nvar-1))
      sum.smc <- sum(smc(r,covar=covar))
      lambda.6 <- (Vt +sum.smc-sum(diag(r)))/Vt
    omega_h <- gsq/Vt	
    om.tot <- (Vt-uniq)/Vt
      om.limit <- gsq/(Vt-uniq)     
    alpha <- ((Vt-Vitem)/Vt)*(nvar/(nvar-1))
    
    Loads<- cbind(gload,group)
    R2<- fsi(f=Loads,phi=NULL,r=r)
    
 #find the subset omegas   (taken from omega)
     omg <- omgo <- omt<-  rep(NA,nfactors+1)
     sub <- apply(loads,1,function(x) which.max(abs(x[2:(nfactors+1)])))  #find the items
     grs <- 0
     for(group in (1:nfactors)) {
     groupi <- which(sub==group)
     if(length(groupi) > 0) {
      Vgr <- sum(r[groupi,groupi])
      gr <- sum(loads[groupi,(group+1)])  #?
      grs <- grs + gr^2
      omg[group+1] <- gr^2/Vgr
      omgo[group+1] <- sum(loads[groupi,1])^2/Vgr
      omt[group+1] <- (gr^2+ sum(loads[groupi,1])^2)/Vgr
     }
     omgo[1] <- sum(loads[,1])^2/sum(r)  #omega h
     omg[1] <- grs/sum(r)  #omega of subscales
     omt[1] <- om.tot 
     om.group <- data.frame(total=omt,general=omgo,group=omg)
     rownames(om.group) <- colnames(loads)[1:(nfactors+1)]
     om.group <- as.matrix(om.group)

     }   

 result <- list(
 Call=cl,
 omega_h=omega_h,
 omega_lim=om.limit,
 omega_total=om.tot,
 alpha=alpha,
 R2=R2,Phi=Phi,
 omega.group=om.group, 
 loadings=loads,
 communality =h2)	   
class(result) <- c("psych","CFA.bifactor")
return(result)
}

if(FALSE) {
Loads<- cbind(gload,group)
sigma <-  Loads %*% Phi %*% t(Loads)
FDs <- diag(Phi %*% t(lambda) %*% solve(sigma) %*% lambda %*% Phi)^0.5
FD <- diag(PHI %*% t(group) %*% solve(r) %*% group %*% Phi)^0.5
}