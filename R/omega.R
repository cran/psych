"omega" <-
function(m,nfactors=3,pc="pa",...) {
      #m is a correlation matrix
      #nfactors is the number of factors to extract
      require(GPArotation)
      nvar <-dim(m)[2]
      gf<-schmid(m,nfactors,pc,...)
      Vt <- sum(m)   #find the total variance in the scale
      Vitem <-sum(diag(m)) #
      gload <- gf$sl[,1]
      gsq <- (sum(gload))^2
      alpha <- ((Vt-Vitem)/Vt)*(nvar/(nvar-1))
      omega <- list(omega= gsq/Vt,alpha=alpha,schmid=gf)
      }

