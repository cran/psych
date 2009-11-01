 "guttman" <- 
 function(r,key=NULL,digits=2) {
 cl <- match.call()
 nvar <- dim(r)[2]
if(dim(r)[1] != dim(r)[2]) {r <- cor(r,use="pairwise")}  else {
if(!is.matrix(r)) r <- as.matrix(r)
r <- cov2cor(r)}  #make sure it is a correlation matrix not a covariance or data matrix
      if(is.null(colnames(r))) {  rownames(r) <- colnames(r) <- paste("V",1:nvar,sep="") }
      if (!is.null(key)) { key <- as.vector(key)
                          r <- diag(key) %*% r %*% diag(key)
                        
                           flip <- FALSE   #we do this if we specify the key
                           } else {key <- rep(1,nvar) }
      
      m <- (1-r)/2
      diag(m) <- 1
      m.names <- colnames(r)
      colnames(m) <- m.names   #flip items if we choose to do so
       
 

       signkey <- strtrim(key,1)
       signkey[signkey=="1"] <- ""
        m.names <- paste(m.names,signkey,sep="")
        colnames(m) <- rownames(m) <- m.names 
            
 if(nvar < 3) {message("These estimates are not really meaningful if you have less than 3 items, \n Try running the alpha function instead")
            stop()}
 worst <- ICLUST(r,2,plot=FALSE)
 w.keys <- worst$p.sorted$cluster
 
 best <- ICLUST(m,2,plot=FALSE)
 keys <- matrix(rep(0,nvar*2),ncol=2)
 b.keys <- best$p.sorted$cluster
 
 m1 <- r
 diag(m1) <- 0
 best.kmeans <- kmeans(m,2,nstart=10)
 keys.kmean <- matrix(rep(0,nvar*2),ncol=2)
 for(i in 1:nvar) {
 keys.kmean[i,best.kmeans$cluster[i]] <- 1 }  
 
  f1 <- factor.pa(r,SMC=FALSE)  #one factor solution
  load <- f1$loadings
   ord.load <- order(load)
    key.fa <- matrix(rep(0,nvar*2),ncol=2)
    for (i in 1:nvar) {
    key.fa[ord.load[i],1] <- i %% 2
    key.fa[ord.load[i],2] <- 1 -  key.fa[ord.load[i],1] }
  
  
    
  f2 <- factor.pa(r,2,SMC=FALSE)  #two factor solution
    load <- f2$loadings
     key.fa2 <- matrix(rep(0,nvar*2),ncol=2)
   
    key.fa2[,1] <- (load[,1] > load[,2]) + 0
    key.fa2[,2 ] <- 1- key.fa2[,1]
    
ev <-eigen(r)$values
e <- ev[1]
alpha.pc <- (1-1/e) * nvar/(nvar-1)
#alpha.pc2 <- (1-1/ev[2]) * nvar/(nvar-1)

r.pc <-  2*ev[1]/(ev[1]+ev[2])-1
r.pc <- r.pc * alpha.pc #attenuate the correlation
beta.pc <- 2 * r.pc/(1+r.pc)

Vt <-  sum.r <- sum(r)
 tr.r <- tr(r)
 lambda.1 <- 1 - tr.r/Vt
 
 off <- r
diag(off) <- 0
sum.off <- sum(off)
sumsq.off <- sum(off^2)
lambda.2 <- (sum.off+ sqrt(sumsq.off*nvar/(nvar-1)))/Vt
lambda.3 <- nvar * lambda.1/(nvar-1)
sum.smc <- sum(smc(r))
lambda.6 <-(sum.r+sum.smc-sum(diag(r)))/Vt

c.co <- colSums(r^2)-diag(r^2)
c.co.max <- max(c.co)
lambda.5 <- lambda.1 + 2*sqrt(c.co.max)/Vt
lambda.5p <- lambda.1 +(nvar)/(nvar-1)*  2*sqrt(c.co.max)/Vt
 
keys <- cbind(w.keys,b.keys,keys.kmean,key.fa,key.fa2)
try(colnames(keys) <- c("IC1","IC2","ICr1","ICr2","K1","K2","F1","F2","f1","f2"))
 covar <- t(keys) %*% r %*% keys    #matrix algebra is our friend
 var <- diag(covar)
 sd.inv <- 1/sqrt(var)
 ident.sd <- diag(sd.inv,ncol = length(sd.inv))
 cluster.correl <- ident.sd %*% covar  %*% ident.sd
 beta <- abs(cluster.correl[2,1]) *2 /(1+abs(cluster.correl[2,1]))
 glb1 <- cluster.correl[3,4] *2 /(1+cluster.correl[3,4])
 glb2 <- cluster.correl[5,6] * 2/(1+ cluster.correl[5,6] )
 glb3 <- cluster.correl[7,8] * 2/(1+cluster.correl[7,8])
 beta.fa <- cluster.correl[9,10] * 2/(1+cluster.correl[9,10])
 glb.max <- max(glb1,glb2,glb3)
 sum.smc <- sum(smc(r))

 gamma <- (sum.r+sum.smc-sum(diag(r)))/Vt
 tenberg <- tenberge(r,digits=digits)
 result <- list(lambda.1=round(lambda.1,digits),lambda.2=round(lambda.2,digits),lambda.3=round(lambda.3,digits),lambda.4 =round(glb.max,digits),lambda.5 = round(lambda.5,digits),lambda.5p = round(lambda.5p,digits),lambda.6=round(lambda.6,digits),beta = round(beta,digits),beta.factor = round(beta.fa,digits),alpha.pc = round(alpha.pc,digits),
 glb.IC =round(glb1,digits),glb.Km = round(glb2,digits), glb.Fa =round(glb3,digits), keys=keys, tenberge=tenberg,r.pc=r.pc,beta.pc=beta.pc,Call=cl)
 class(result) <- c("psych","guttman")
 return(result)
}
