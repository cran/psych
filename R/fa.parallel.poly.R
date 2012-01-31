 #parallel analysis of polychoric factor analysis
 "fa.parallel.poly" <- 
function(x,n.iter=10,SMC=TRUE,fm="minres") { 
p <- .05
 cl <- match.call()

n.obs <- dim(x)[1]
	tx <- table(as.matrix(x))
	nx <- dim(tx)[1] 
	if(dim(tx)[1] ==2) {tet <- tetrachoric(x)
	    typ = "tet"} else {tet <- polychoric(x)
	    typ = "poly"}
 r <- tet$rho
 f <- fa(r,n.obs=n.obs,SMC = SMC,covar=FALSE,fm=fm) #call fa with the appropriate parameters
 f$Call <- cl
 fl <- f$loadings  #this is the original
 nvar <- dim(fl)[1]

 e.values <- list(pc =list(),fa = list())
 replicates <- list()
 rep.rots <- list()
 for (trials in 1:n.iter) {
 xs <- matrix(sample(nx,n.obs*nvar,replace=TRUE),ncol=nvar)
  if(typ=="poly") {tets <- polychoric(xs)} else {tets <- tetrachoric(xs)}
  r <- tets$rho
  values.samp <- eigen(tets$rho)$values
  e.values[["pc"]][[trials]] <- values.samp
   					
   				
  fs <- fa(r,n.obs=n.obs,SMC = SMC,covar=FALSE,fm=fm) #call fa with the appropriate parameters
     e.values[["fa"]][[trials]] <- fs$values
     
} 


ei.pc <-describe(matrix(unlist(e.values$pc),ncol=nvar,byrow=TRUE))  #eigen values of pcs
ei.fa <- describe(matrix(unlist(e.values$fa),ncol=nvar,byrow=TRUE)) #eigen values of fa
fa.test <- which(!(f$values > ei.fa$mean))[1]-1
pc.test <- which(!(f$e.values > ei.pc$mean))[1] -1
#e.stats <- list(fa.values=f$values,pc.values=f$e.values,pc.sim=ei.pc,fa=ei.fa,nfact=fa.test,ncomp = pc.test)


#results <- list(fa = f,rho=tet$rho,tau=tet$tau,n.obs=n.obs,Call= cl,e.values=e.values,e.stats=e.stats)
results <- list(rho=tet$rho,tau=tet$tau,n.obs=n.obs,Call= cl,fa.values=f$values,pc.values=f$e.values,pc.sim=ei.pc,fa.sim=ei.fa,nfact=fa.test,ncomp = pc.test)
class(results) <- c("psych","parallel")
plot.poly.parallel(results)
return(results)
 }
 #written May 8 2011
 #corrected December 27, 2011 to pass the fm correctly (had always forced minres)  
 
 
 "plot.poly.parallel" <-
 function(x,show.legend=TRUE,...) {
 e.values <- x$pc.values
 values <-  x$fa.values
pcm <- x$pc.sim$mean
fam <- x$fa.sim$mean

 plot(e.values,type="b", main = "Eigen values of tetrachoric/polychoric matrix",ylab="Eigen values of original  and simulated factors and components",ylim=c(0,max(e.values)) ,xlab="Factor Number",pch=4,col="blue")
 points(values,type ="b",pch=2,col="blue")
 points(pcm,type="l",lty="dotted",pch=2,col="red")
 points(fam,type="l",lty="dashed",pch=2,col="red")
 if(show.legend) {

       legend("topright", c("PC  Actual Data", " PC  Simulated Data","FA  Actual Data", " FA  Simulated Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
       text.col = "green4", lty = c("solid","dotted","solid","dotted"),
       merge = TRUE, bg = 'gray90')
   }
   
fa.test <- which(!(values > fam))[1]-1
pc.test <- which(!(e.values > pcm))[1] -1
cat("Parallel analysis suggests that ")
cat("the number of factors = ",fa.test, " and the number of components = ",pc.test,"\n")
 }