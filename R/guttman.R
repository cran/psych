#revised December 2, 2013 to take advantage of the splitHalf function
 "guttman" <- 
 function(r,key=NULL) {
 cl <- match.call()
 .Deprecated("splitHalf",msg="Guttman has been deprecated.  The use of the splitHalf function is recommended") 
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
# beta <- ICLUST(r,1,plot=FALSE)$beta
# worst <- ICLUST(r,2,plot=FALSE)
# w.keys <- worst$p.sorted$cluster
 
 #the following was a crude attempt at finding the best
 #this has been replaced with calling splitHalf
 best <- splitHalf(r)
 # best <- ICLUST(m,2,plot=FALSE,SMC=FALSE)
 #best <- ICLUST(m,2,plot=FALSE)
 #keys <- matrix(rep(0,nvar*2),ncol=2)
 #b.keys <- best$p.sorted$cluster
 
# m1 <- r
 #diag(m1) <- 0
# best.kmeans <- kmeans(m,2,nstart=10)
#keys.kmean <- matrix(rep(0,nvar*2),ncol=2)
# for(i in 1:nvar) {
# keys.kmean[i,best.kmeans$cluster[i]] <- 1 }  
 
  f1 <- fa(r,SMC=FALSE)  #one factor solution
  load <- f1$loadings
   ord.load <- order(load)
    key.fa <- matrix(rep(0,nvar*2),ncol=2)
    for (i in 1:nvar) {
    key.fa[ord.load[i],1] <- i %% 2
    key.fa[ord.load[i],2] <- 1 -  key.fa[ord.load[i],1] }
  
  
    
  f2 <- fa(r,2,SMC=FALSE)  #two factor solution
    load <- f2$loadings
     key.fa2 <- matrix(rep(0,nvar*2),ncol=2)
   
    key.fa2[,1] <- (abs(load[,1]) > abs(load[,2])) + 0
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
#
#this next section is a complete kludge meant to find the most similar splits
#a better way is to use the glb function of Andreas Moeltner
#revised February 11 to implement equation 51 of Guttman, not 51'
#
#all of this has been deleted as of December, 2013 to just us splitHalf
#keys <- cbind(w.keys,b.keys,keys.kmean,key.fa,key.fa2)
#try(colnames(keys) <- c("IC1","IC2","ICr1","ICr2","K1","K2","F1","F2","f1","f2"))
 #covar <- t(keys) %*% r %*% keys    #matrix algebra is our friend
# var <- diag(covar)
# sd.inv <- 1/sqrt(var)
# ident.sd <- diag(sd.inv,ncol = length(sd.inv))
# cluster.correl <- ident.sd %*% covar  %*% ident.sd
 #beta <- abs(cluster.correl[2,1]) *2 /(1+abs(cluster.correl[2,1])) #worst split 
# beta <- 2 * (1-2/(2+abs(cluster.correl[2,1])))
 #glb1 <- cluster.correl[3,4] *2 /(1+cluster.correl[3,4])
# glb2 <- cluster.correl[5,6] * 2/(1+ cluster.correl[5,6] )
# glb3 <- cluster.correl[7,8] * 2/(1+cluster.correl[7,8])
#Vtcl1 <- covar[3,3]+ covar[4,4] + 2 * covar[3,4]
#Vtcl2 <- covar[5,5]+ covar[6,6] + 2 * covar[5,6]
#Vtcl3 <- covar[7,7]+covar[8,8] + 2 * covar[7,8]
#glbIC <- 2*(1-(covar[3,3]+ covar[4,4])/Vtcl1  )
#glb2 <- 2*(1-(covar[5,5]+ covar[6,6])/Vtcl2  )
#glb3 <- 2*(1-(covar[7,7]+ covar[8,8])/Vtcl3  )
#beta.fa <- cluster.correl[9,10] * 2/(1+cluster.correl[9,10])
# glb.max <- max(glbIC,glb2,glb3)
 sum.smc <- sum(smc(r))
 glb <- glb.fa(r)$glb
beta <- best$minrb
if(beta < 0) beta <- 0
 gamma <- (sum.r+sum.smc-sum(diag(r)))/Vt
 tenberg <- tenberge(r)
 result <- list(lambda.1=lambda.1,lambda.2=lambda.2,lambda.3=lambda.3,lambda.4 =best$maxrb,lambda.5 = lambda.5,lambda.5p = lambda.5p,lambda.6=lambda.6,alpha.pc = alpha.pc,
 glb=glb, tenberge=tenberg,r.pc=r.pc,beta.pc=beta.pc,beta=beta,Call=cl)
 class(result) <- c("psych","guttman")
 return(result)
}
