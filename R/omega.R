"omegah" <-
function(m,nfactors=3,fm="minres",key=NULL,flip=TRUE, digits=2,title="Omega",sl=TRUE,labels=NULL, plot=TRUE,n.obs=NA,rotate="oblimin",Phi = NULL,option="equal",covar=FALSE,...) {
      #m is a correlation matrix, or if not, the correlation matrix is found
      #nfactors is the number of factors to extract
      #key allows items to be reversed scored  if desired
      #if Phi is not null, this implies that we have been given a factor matrix  -- added May 30, 2010
     if(!requireNamespace('GPArotation') && (rotate !="cluster")) {stop("I am sorry, you need to have the  GPArotation package installed")}
      cl <- match.call() 
      nvar <- dim(m)[2]
      raw.data <- NULL
      if(is.null(Phi)) {   #the normal case is to do the factor analysis of the raw data or the correlation matrix
      if(dim(m)[1] != dim(m)[2]) {
                            n.obs <- dim(m)[1]
                            m <- as.matrix(m)
                            raw.data <- m #added 9/1/14
                            if(covar) {m <- cov(m,use="pairwise") } else {m <- cor(m,use="pairwise")}
                            } else {
                            if(!covar) m <- cov2cor(as.matrix(m))    #make sure it is a correlation matrix not a covariance or data matrix (if we change this, we will need to change the calculation for omega later)
                           }
                            
      if(is.null(colnames(m))) {  rownames(m) <- colnames(m) <- paste("V",1:nvar,sep="") }
       m.names <- colnames(m)
      
      if (!is.null(key)) { m <- diag(key) %*% m %*% diag(key)
                           colnames(m) <- m.names   #flip items if we choose to do so
                           flip <- FALSE   #we do this if we specify the key
                           } else {key <- rep(1,nvar) }
             signkey <- strtrim(key,1)
             signkey[signkey=="1"] <- ""
             m.names <- paste(m.names,signkey,sep="")
             colnames(m) <- rownames(m) <- m.names
      if ((nvar < 6) && (fm =="mle") ) {message(paste("In omega, 3 factors are too many for ",nvar," variables using mle.  Using minres instead",sep=""))
       fm <- "minres"} 
       } else { m.names <- rownames(m) }   #add the names if we have a factor input 
       
      gf <-schmid(m,nfactors,fm,digits,rotate=rotate,n.obs=n.obs,Phi=Phi,option=option,covar=covar, ...)
      if(!is.null(Phi)) { model <- m
                         nfactors <- dim(model)[2]
                          m <- factor.model(model,Phi=Phi,U2=FALSE)
                          nvar <- dim(m)[2]
                          if(is.null(rownames(m))) {colnames(m) <- rownames(m) <- paste("V",1:nvar)}
                          }
                         
      gload <- gf$sl[,1]
      if (flip) {       #should we think about flipping items ?
       			key <- sign(gload)
      		 	key[key==0] <- 1     # a rare and weird case where the gloading is 0 and thus needs not be flipped
      			 if (sum(key) < nvar) {  #some items have negative g loadings and should be flipped  
            		 m <- diag(key) %*% m %*% diag(key)  #this is just flipping the correlation matrix so we can calculate alpha
            		 gf$sl[,1:(nfactors+1)] <-  diag(key) %*% gf$sl[,1:(nfactors+1)] 
            		 signkey <- strtrim(key,1)
            		 signkey[signkey=="1"] <- ""
            		 m.names <- paste(m.names,signkey,sep="")
            		 colnames(m) <- rownames(m) <- m.names
            		 rownames(gf$sl) <- m.names
          				
     		 
            }
          }
          Vt <- sum(m)   #find the total variance in the scale
     	  Vitem <- sum(diag(m)) 
     	  gload <- gf$sl[,1]	
       
     
      gsq <- (sum(gload))^2
      uniq <- sum(gf$sl[,"u2"])
      
      if((nfactors == 1) && (fm=="pc")) {gsq <- Vt - uniq
                                        warning("omega_h is not meaningful for a principal components analysis with one component")}   #weird condition when using fm=pc and 1 factor
      om.tot <- (Vt-uniq)/Vt
      om.limit <- gsq/(Vt-uniq)  
      alpha <- ((Vt-Vitem)/Vt)*(nvar/(nvar-1))
      sum.smc <- sum(smc(m,covar=covar))
      lambda.6 <- (Vt +sum.smc-sum(diag(m)))/Vt
      if (!is.null(digits)) {omega <-list(omega_h= gsq/Vt,alpha=alpha,lambda.6 = lambda.6,omega.tot =om.tot,schmid=gf ,key = key,title=title)
            dg <-max(digits-1,1)} else {
              omega <- list(omega_h= gsq/Vt,alpha=alpha,omega.tot=om.tot,schmid=gf,key=key,title=title)
      dg <- 1}
      ev <- colSums(gf$sl[,1:(nfactors+1)]^2)
      ECV <- ev[1]/sum(ev)
      omega.stats <- factor.stats(m,gf$sl[,1:(nfactors+1)],n.obs=n.obs)
      general.stats <- factor.stats(m,as.matrix(gf$sl[,1]),n.obs=n.obs)  #just get fit for the general factor
     if (nfactors<2) plot <- FALSE
    # if(require(Rgraphviz) && plot) {omega.model <-omega.graph(omega,title=title,sl=sl,labels=labels,digits=dg) } else {omega.model <- omega.sem(omega,sl=sl)}
   
     omega.model <- omega.sem(omega,sl=sl)
     
     #find the subset omegas
     omg <- omgo <- omt<-  rep(NA,nfactors+1)
     sub <- apply(gf$sl,1,function(x) which.max(abs(x[2:(nfactors+1)])))
     grs <- 0
     for(group in( 1:nfactors)) {
     groupi <- which(sub==group)
     if(length(groupi) > 0) {
      Vgr <- sum(m[groupi,groupi])
      gr <- sum(gf$sl[groupi,(group+1)])
      grs <- grs + gr^2
      omg[group+1] <- gr^2/Vgr
      omgo[group+1] <- sum(gf$sl[groupi,1])^2/Vgr
      omt[group+1] <- (gr^2+ sum(gf$sl[groupi,1])^2)/Vgr
     }
     omgo[1] <- sum(gf$sl[,1])^2/sum(m)  #omega h
     omg[1] <- grs/sum(m)  #omega of subscales
     omt[1] <- om.tot 
     om.group <- data.frame(total=omt,general=omgo,group=omg)
     rownames(om.group) <- colnames(gf$sl)[1:(nfactors+1)]

     if(!is.null(raw.data)) {scores <- raw.data %*%  omega.stats$weights} else {scores<- NULL} 
     }
     omega <- list(omega_h= gsq/Vt,omega.lim = om.limit,alpha=alpha,omega.tot=om.tot,G6=lambda.6,schmid=gf,key=key,stats = omega.stats,ECV=ECV,gstats = general.stats,call=cl,title=title,R = m,model=omega.model,omega.group=om.group,scores=scores)

      class(omega) <- c("psych","omega")
     if(plot)  omega.diagram(omega,main=title,sl=sl,labels=labels,digits=dg)
      return(omega)
      }
#April 4, 2011  added a check for fm=pc and nfactors == 1 to solve problem of omega_h < omega_t -- probably not a good idea.  removed
#January 9, 2014  added omega scores if the raw data are given


"omega" <- 
function(m,nfactors=3,fm="minres",n.iter=1,p=.05,poly=FALSE,key=NULL,flip=TRUE, digits=2,title="Omega",sl=TRUE,labels=NULL, plot=TRUE,n.obs=NA,rotate="oblimin",Phi = NULL,option="equal",covar=FALSE,...) {
 cl <- match.call()
  if(is.data.frame(m) || is.matrix(m)) {if(dim(m)[1] == dim(m)[2] ) {if(is.na(n.obs) && (n.iter>1)) stop("You must specify the number of subjects if giving a correlation matrix")
                                # if(!require(MASS)) stop("You must have MASS installed to simulate data from a correlation matrix")
                                 }
                  }
 if(!is.data.frame(m) && !is.matrix(m)) {  n.obs=m$n.obs
                if(poly) {
                      pol <- list(rho=m$rho,tau = m$tau,n.obs=m$n.obs)
                        
                           m <- m$rho
                          } else { m <- m$R 
                          }
                          } else { #new data

          if(poly) { pol <- polychoric(m)
                     m <- pol$rho
                     n.obs <- pol$n.obs}  }       
            om <- omegah(m=m,nfactors=nfactors,fm=fm,key=key,flip=flip, digits=digits,title=title,sl=sl,labels=labels, plot=plot,n.obs=n.obs,rotate=rotate,Phi = Phi,option=option,covar=covar,...) #call omega with the appropriate parameters
    

       
 if(is.na(n.obs) ) {n.obs <- om$stats$n.obs} 
 replicates <- list()
 if(n.iter > 1) {for (trials in 1:n.iter) {
 if(dim(m)[1] == dim(m)[2]) {#create data sampled from multivariate normal with correlation
                                      nvar <- dim(m)[1]
                                      #mu <- rep(0, nvar)
                                     # m <- mvrnorm(n = n.obs, mu, Sigma = m, tol = 1e-06, empirical = FALSE)
 #the next 3 lines replaces mvrnorm (taken from mvrnorm, but without the checks)
                                      eX <- eigen(m)
                                      m <- matrix(rnorm(nvar * n.obs),n.obs)
                                      m <-  t(eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(m))
                            } else {m <- m[sample(n.obs,n.obs,replace=TRUE),]}
    if(poly) {pol <- polychoric(m)
            oms <- omegah(m=pol$rho,nfactors=nfactors,fm=fm,key=key,flip=flip, digits=digits,title=title,sl=sl,labels=labels, plot=plot,n.obs=pol$n.obs,rotate=rotate,Phi = Phi,option=option,...) #call omega with the appropriate parameters
            } else {
 oms <- omegah(m=m,nfactors=nfactors,fm=fm,key=key,flip=flip, digits=digits,title=title,sl=sl,labels=labels, plot=plot,n.obs=n.obs,rotate=rotate,Phi = Phi,option=option,...) #call omega with the appropriate parameters
 }                         
 # oms <-omegah(m=m,nfactors=nfactors,fm=fm,key=key,flip=flip, digits=digits,title=title,sl=sl,labels=labels, plot=plot,n.obs=n.obs,rotate=rotate,Phi = Phi,option=option,...) #call fa with the appropriate parameters
  replicates[[trials]] <- list(omega=oms$omega_h,alpha=oms$alpha,omega.tot=oms$omega.tot,G6=oms$G6,omega.lim=oms$omega.lim)
  }
  
replicates <- matrix(unlist(replicates),ncol=5,byrow=TRUE)
z.replicates <- cbind(fisherz(replicates[,1:4]),replicates[,5])  #convert to z scores
means <- colMeans(z.replicates,na.rm=TRUE)
sds <-  apply(z.replicates,2,sd,na.rm=TRUE)
ci.lower <-  means + qnorm(p/2) * sds
ci.upper <- means + qnorm(1-p/2) * sds
ci <- data.frame(lower = ci.lower,upper=ci.upper)
ci <- rbind(fisherz2r(ci[1:4,]),ci[5,])
rownames(ci) <- c("omega_h","alpha","omega_tot","G6","omega_lim")
colnames(replicates) <- names(means) <- names(sds) <- rownames(ci) 
conf <- list(means = means,sds = sds,ci = ci,Call= cl,replicates=replicates)
om$Call=cl
results <- list(om = om,ci=conf) } else {om$call=cl
                                         if(poly) {om$rho <- pol$rho
                                                   om$tau <- pol$tau
                                                   om$n.obs <- pol$n.obs
                                                   }
                                         results <- om}
class(results) <- c("psych","omega")
return(results)
 }
 #written April 25, 2011
 #adapted May 12, 2011 to be the primary version of omega
 


