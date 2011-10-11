"omega.bifactor" <- function(r,nfactors=3,rotate="bifactor",n.obs = NA,flip=TRUE,key=NULL,title="Omega from bifactor",...) {
 cl <- match.call()
 if(dim(r)[2] != dim(r)[1]) {n.obs <- dim(r)[1]
                             r <- cor(r,use="pairwise")}
 nvar <- dim(r)[2]
  if(is.null(colnames(r))) {  rownames(r) <- colnames(r) <- paste("V",1:nvar,sep="") }
       r.names <- colnames(r)
 
 if (!is.null(key)) { r <- diag(key) %*% r %*% diag(key)
                           colnames(r) <- r.names   #flip items if we choose to do so
                           flip <- FALSE   #we do this if we specify the key
                           } else {key <- rep(1,nvar) }
                           
 f <- fa(r,nfactors=nfactors,rotate=rotate,n.obs = n.obs) 
 if (flip) {       #should we think about flipping items ?
       			key <- sign(f$loadings[,1])
      		 	key[key==0] <- 1     # a rare and weird case where the gloading is 0 and thus needs not be flipped
      			 if (sum(key) < nvar) {  #some items have negative g loadings and should be flipped  
            		 r <- diag(key) %*% r %*% diag(key)  #this is just flipping the correlation matrix so we can calculate alpha
            		f$loadings <-  diag(key) %*% f$loadings
            		 signkey <- strtrim(key,1)
            		 signkey[signkey=="1"] <- ""
            		 r.names <- paste(r.names,signkey,sep="")
            		 colnames(r) <- rownames(r) <- r.names
            		 rownames(f$loadings) <- r.names
          				
     		 
            }
          }
          Vt <- sum(r)   #find the total variance in the scale
     	  Vitem <- sum(diag(r)) 
     	  gload <- f$loadings[,1]	
 
 
       gsq <- (sum(gload))^2
      uniq <- nvar - tr(f$loadings %*% t(f$loadings))
       
      om.tot <- (Vt-uniq)/Vt
      om.limit <- gsq/(Vt-uniq)  
      alpha <- ((Vt-Vitem)/Vt)*(nvar/(nvar-1))
      sum.smc <- sum(smc(r))
      lambda.6 <- (Vt +sum.smc-sum(diag(r)))/Vt
      omega_h= gsq/Vt
      omega <-list(omega_h= omega_h,alpha=alpha,G6 = lambda.6,omega.tot =om.tot ,key = key,title=title,f=f)
      class(omega) <- c("psych","bifactor")
      return(omega)
      }
      
      
      

