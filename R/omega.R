"omega" <-
function(m,nfactors=3,fm="minres",key=NULL,flip=TRUE, digits=2,title="Omega",sl=TRUE,labels=NULL, plot=TRUE,n.obs=NA,rotate="oblimin",...) {
      #m is a correlation matrix, or if not, the correlation matrix is found
      #nfactors is the number of factors to extract
      #key allows items to be reversed scored  if desired
     if(!require(GPArotation) && (rotate !="cluster")) {stop("I am sorry, you need to have the  GPArotation package installed")}
      cl <- match.call()
      nvar <- dim(m)[2]
      if(dim(m)[1] != dim(m)[2]) {
                            n.obs <- dim(m)[1]
                            m <- cor(m,use="pairwise")} else {
                            m <- cov2cor(as.matrix(m))    #make sure it is a correlation matrix not a covariance or data matrix (if we change this, we will need to change the calculation for omega later)
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
       gf<-schmid(m,nfactors,fm,digits,rotate=rotate,n.obs=n.obs, ...)

      Vt <- sum(m)   #find the total variance in the scale
      Vitem <-sum(diag(m)) #
      gload <- gf$sl[,1]
      if (flip) {       #should we think about flipping items ?
        key <- sign(gload)
        key[key==0] <- 1     # a rare and weird case where the gloading is 0 and thus needs not be flipped
        if (sum(key) < nvar) {  #some items have negative g loadings and should be flipped  
            
             if(dim(m)[1] != dim(m)[2]) m <- cor(m,use="pairwise")
             m <- diag(key) %*% m %*% diag(key)
           
             signkey <- strtrim(key,1)
             signkey[signkey=="1"] <- ""
             m.names <- paste(m.names,signkey,sep="")
             colnames(m) <- rownames(m) <- m.names
             
             gf<-schmid(m,nfactors,fm,digits=digits,n.obs=n.obs,...)
      Vt <- sum(m)   #find the total variance in the scale
      Vitem <-sum(diag(m)) #
      gload <- gf$sl[,1]
            }
          }
      gsq <- (sum(gload))^2
      uniq <- sum(gf$sl[,(nfactors+3)])
      om.tot <- (Vt-uniq)/Vt
      om.limit <- gsq/(Vt-uniq)  
      alpha <- ((Vt-Vitem)/Vt)*(nvar/(nvar-1))
      sum.smc <- sum(smc(m))
      lambda.6 <-(Vt +sum.smc-sum(diag(m)))/Vt
      if (!is.null(digits)) {omega <-list(omega_h= gsq/Vt,alpha=alpha,lambda.6 = lambda.6,omega.tot =om.tot,schmid=gf ,key = key,title=title)
      dg <-max(digits-1,1)} else {
      omega <- list(omega_h= gsq/Vt,alpha=alpha,omega.tot=om.tot,schmid=gf,key=key,title=title)
      dg <- 1}
      omega.stats <- factor.stats(m,gf$sl[,1:(nfactors+1)])
     if (nfactors<2) plot <- FALSE
    # if(require(Rgraphviz) && plot) {omega.model <-omega.graph(omega,title=title,sl=sl,labels=labels,digits=dg) } else {omega.model <- omega.sem(omega,sl=sl)}
    omega.model <- omega.sem(omega,sl=sl)
     
     omega <- list(omega_h= gsq/Vt,omega.lim = om.limit,alpha=alpha,omega.tot=om.tot,G6=lambda.6,schmid=gf,key=key,stats = omega.stats,call=cl,title=title,model=omega.model)

      class(omega) <- c("psych","omega")
      omega.diagram(omega,main=title,sl=sl,labels=labels,digits=dg)
      return(omega)
      }

