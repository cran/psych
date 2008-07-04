"omega" <-
function(m,nfactors=3,pc="mle",key=NULL,flip=TRUE, digits=2,title="Omega",sl=TRUE,labels=NULL, plot=TRUE,rotate="oblimin",...) {
      #m is a correlation matrix, or if not, the correlation matrix is found
      #nfactors is the number of factors to extract
      #key allows items to be reversed scored  if desired
      if(!require(GPArotation)) {stop("I am sorry, you need to have the  GPArotation package installed")}
      nvar <-dim(m)[2]
      if(dim(m)[1] != dim(m)[2]) m <- cor(m,use="pairwise")
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
      if ((nvar < 6) && (pc=="mle") ) {warning(paste("3 factors is too many for ",nvar," variables using mle.  Using pa instead",sep=""))
       pc <- "pa"} 
       gf<-schmid(m,nfactors,pc,digits,rotate=rotate, ...)

      Vt <- sum(m)   #find the total variance in the scale
      Vitem <-sum(diag(m)) #
      gload <- gf$sl[,1]
      if (flip) {       #should we think about flipping items ?
        key <- sign(gload)
        if (sum(key) < nvar) {  #sum items have negative g loadings and should be flipped  
            
             if(dim(m)[1] != dim(m)[2]) m <- cor(m,use="pairwise")
             m <- diag(key) %*% m %*% diag(key)
           
             signkey <- strtrim(key,1)
             signkey[signkey=="1"] <- ""
             m.names <- paste(m.names,signkey,sep="")
             colnames(m) <- rownames(m) <- m.names
             
             gf<-schmid(m,nfactors,pc,digits=digits,...)
      Vt <- sum(m)   #find the total variance in the scale
      Vitem <-sum(diag(m)) #
      gload <- gf$sl[,1]
            }
          }
      gsq <- (sum(gload))^2
      alpha <- ((Vt-Vitem)/Vt)*(nvar/(nvar-1))
      if (!is.null(digits)) {omega <-list(omega_h= round(gsq/Vt,digits),alpha=round(alpha,digits),schmid=gf ,key = key) } else {
      omega <- list(omega_h= gsq/Vt,alpha=alpha,schmid=gf,key=key) }
      if (plot) {omega.graph(omega,title=title,sl=sl,labels=labels) }
      return(omega)
      }

