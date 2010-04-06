# an attempt at finding the worst and best splits, beta is worst split (from ICLUST)

"glb" <-
function(r,key=NULL) {
 nvar <- dim(r)[2]
if(dim(r)[1] != dim(r)[2]) {r <- cor(r,use="pairwise")} else {r <- cov2cor(r)}  #make sure it is a correlation matrix not a covariance or data matrix
      if(is.null(colnames(r))) {  rownames(r) <- colnames(r) <- paste("V",1:nvar,sep="") }
      m <- (1-r)/2
      diag(m) <- 1
       m.names <- colnames(m)
       
 
 if (!is.null(key)) { m <- diag(key) %*% m %*% diag(key)
                           colnames(m) <- m.names   #flip items if we choose to do so
                           flip <- FALSE   #we do this if we specify the key
                           } else {key <- rep(1,nvar) }
       signkey <- strtrim(key,1)
             signkey[signkey=="1"] <- ""
             m.names <- paste(m.names,signkey,sep="")
             colnames(m) <- rownames(m) <- m.names             
 worst <- ICLUST(r,2,plot=FALSE)
  keys <- worst$p.sorted$cluster
 
 best <- ICLUST(m,2,plot=FALSE,SMC=FALSE)
 keys <- matrix(rep(0,nvar*2),ncol=2)
 keys <- best$p.sorted$cluster
 
 m1 <- r
 diag(m1) <- 0
 best.kmeans <- kmeans(m,2,nstart=10)
 keys.kmean <- matrix(rep(0,nvar*2),ncol=2)
 for(i in 1:nvar) {
 keys.kmean[i,best.kmeans$cluster[i]] <- 1 }  
 
  f1 <- fa(r)  #one factor solution
  load <- f1$loadings
   ord.load <- order(load)
    key.fa <- matrix(rep(0,nvar*2),ncol=2)
    for (i in 1:nvar) {
    key.fa[ord.load[i],1] <- i %% 2
    key.fa[ord.load[i],2] <- 1 -  key.fa[ord.load[i],1] }
  
  
    
  f2 <- fa(r,2,SMC=FALSE)  #two factor solution
    load <- f2$loadings
     key.fa2 <- matrix(rep(0,nvar*2),ncol=2)
   
    key.fa2[,1] <- (load[,1] > load[,2]) + 0
    key.fa2[,2 ] <- 1- key.fa2[,1]
    

e <- eigen(r)$values[1]
alpha.pc <- 1-1/e
    
keys <- cbind(worst$p.sorted$cluster,keys,keys.kmean,key.fa,key.fa2)
colnames(keys) <- c("IC1","IC2","ICr1","ICr2","K1","K2","F1","F2","f1","f2")
 covar <- t(keys) %*% r %*% keys    #matrix algebra is our friend
 var <- diag(covar)
 sd.inv <- 1/sqrt(var)
 ident.sd <- diag(sd.inv,ncol = length(sd.inv))
 cluster.correl <- ident.sd %*% covar  %*% ident.sd
 beta <- cluster.correl[2,1] *2 /(1+cluster.correl[2,1])
 glbIC <- cluster.correl[3,4] *2 /(1+cluster.correl[3,4])
 glb2 <- cluster.correl[5,6] * 2/(1+ cluster.correl[5,6] )
 glb3 <- cluster.correl[7,8] * 2/(1+cluster.correl[7,8])
 beta.fa <- cluster.correl[9,10] * 2/(1+cluster.correl[9,10])
 glb.max <- max(glbIC,glb2,glb3)
 sum.smc <- sum(smc(r))
 sum.r <- sum(r)
 gamma <- (sum.r+sum.smc-sum(diag(r)))/sum.r
 tenberg <- tenberge(r)
 result <- list(beta = beta,beta.factor = beta.fa,alpha.pc = alpha.pc, glb.max = glb.max, glb.IC =glbIC,glb.Km = glb2, glb.Fa =glb3, r.smc = gamma,tenberge=tenberg, keys=keys)
 return(result)
}

