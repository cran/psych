"omega" <-
function(m,nfactors=3,pc="mle",key=NULL,...) {
      #m is a correlation matrix, or if not, the correlation matrix is found
      #nfactors is the number of factors to extract
      #key allows items to be reversed scored  if desired
      if(!require(GPArotation)) {stop("I am sorry, you need to have the  GPArotation package installed")}
      nvar <-dim(m)[2]
      if(dim(m)[1] != dim(m)[2]) m <- cor(m,use="pairwise")
      if (!is.null(key)) { m <- diag(key) %*% m %*% diag(key)}   #flip items if we choose to do so}
      if ((nvar < 6) && (pc=="mle") ) {warning(paste("3 factors is too many for ",nvar," variables using mle.  Using pa instead",sep=""))
       pc <- "pa"}
       gf<-schmid(m,nfactors,pc,...)
      Vt <- sum(m)   #find the total variance in the scale
      Vitem <-sum(diag(m)) #
      gload <- gf$sl[,1]
      gsq <- (sum(gload))^2
      alpha <- ((Vt-Vitem)/Vt)*(nvar/(nvar-1))
      omega <- list(omega= gsq/Vt,alpha=alpha,schmid=gf)
      }

