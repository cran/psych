 #parallel analysis of polychoric factor analysis
 "fa.parallel.poly" <- 
function(x,n.iter=10,SMC=TRUE,fm="minres",correct=TRUE,sim=FALSE) { 
p <- .05
 cl <- match.call()

n.obs <- dim(x)[1]
	tx <- table(as.matrix(x))
	nx <- dim(tx)[1] 
	if(dim(tx)[1] ==2) {tet <- tetrachoric(x,correct=correct)
	    typ = "tet"
	    if(sim) { tx.item <- matrix(apply(x,2,table),ncol=ncol(x))
	    px <- matrix(tx.item/colSums(tx.item),ncol=2,byrow=TRUE)
	    } } else {tet <- mixed.cor(x)
	    typ = "poly"}
	    cat("\n")  #to try to clear the progress bar
        flush(stdout())
 r <- tet$rho
 f <- fa(r,n.obs=n.obs,SMC = SMC,covar=FALSE,fm=fm) #call fa with the appropriate parameters
 f$Call <- cl
 fl <- f$loadings  #this is the original
 nvar <- dim(fl)[1]

 e.values <- list(pc =list(),fa = list())
 replicates <- list()
 rep.rots <- list()
 for (trials in 1:n.iter) {
progressBar(trials,n.iter,"fa.parallel.poly")
 #xs <- matrix(sample(nx,n.obs*nvar,replace=TRUE),ncol=nvar) #this does not replicate the p values
 xs <- matrix(apply(x,2,function(y) sample(y,n.obs,replace=TRUE)),ncol=nvar) #do it column wise
# if(typ=="poly") {tets <- polychoric(xs,progress=FALSE)} else {tets <- tetrachoric(xs,correct=correct)}
tets <- polychoric(xs,progress=FALSE)
# tets <- mixed.cor(xs)
  r <- tets$rho
  values.samp <- eigen(tets$rho)$values
  e.values[["pc"]][[trials]] <- values.samp
   					
   				
  fs <- fa(r,n.obs=n.obs,SMC = SMC,covar=FALSE,fm=fm) #call fa with the appropriate parameters
     e.values[["fa"]][[trials]] <- fs$values
     
if(sim && (typ=="tet")) {xsim <- matrix(apply(px,1,function(x) sample(2,n.obs,TRUE,prob=x)),ncol=nvar)
    rho.sim <- tetrachoric(xsim,correct=correct)
    values.sim <- eigen(rho.sim$rho)$values
    e.values[["pc.sim"]][[trials]] <- values.sim
    fsim <- fa(rho.sim$rho,n.obs=n.obs,SMC = SMC,covar=FALSE,fm=fm) #call fa with the appropriate parameters
     e.values[["fa.sim"]][[trials]] <- fsim$values

}
     
} 

cat("\n")  #to try to clear the progress bar

ei.pc <-describe(matrix(unlist(e.values$pc),ncol=nvar,byrow=TRUE))  #eigen values of pcs
ei.fa <- describe(matrix(unlist(e.values$fa),ncol=nvar,byrow=TRUE)) #eigen values of fa
fa.test <- which(!(f$values > ei.fa$mean))[1]-1
pc.test <- which(!(f$e.values > ei.pc$mean))[1] -1
#e.stats <- list(fa.values=f$values,pc.values=f$e.values,pc.sim=ei.pc,fa=ei.fa,nfact=fa.test,ncomp = pc.test)
if(sim && (typ=="tet")) {eis.pc <- describe(matrix(unlist(e.values$pc.sim),ncol=nvar,byrow=TRUE))  #eigen values of pcs                
       eis.fa <- describe(matrix(unlist(e.values$fa.sim),ncol=nvar,byrow=TRUE)) #eigen values of fa
} else {eis.pc <- NULL
        eis.fa <- NULL
}

#results <- list(fa = f,rho=tet$rho,tau=tet$tau,n.obs=n.obs,Call= cl,e.values=e.values,e.stats=e.stats)
results <- list(rho=tet$rho,tau=tet$tau,n.obs=n.obs,Call= cl,fa.values=f$values,pc.values=f$e.values,pc.sim=ei.pc,fa.sim=ei.fa,pcs.sim=eis.pc,fas.sim=eis.fa,nfact=fa.test,ncomp = pc.test)
class(results) <- c("psych","parallel")
cat('\n See the graphic output for a description of the results\n')
plot.poly.parallel(results)
return(results)
 }
 #written May 8 2011
 #corrected December 27, 2011 to pass the fm correctly (had always forced minres)  
 #modified Sept 16, 2013 to use mixed.cor instead of polychoric  (more robust) 
 #modified Oct 2, 2013 to add the ability to do random data as well 
 
 
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
 if(!is.null(x$pcs.sim)) points(x$pcs.sim$mean,type="l",lty="dashed",pch=2,col="red")
 if(!is.null(x$fas.sim)) points(x$fas.sim$mean,type="l",lty="dashed",pch=2,col="red")
 if(show.legend && (!is.null(x$pcs.sim))) {

       legend("topright", c("PC  Actual Data", " PC  Resampled Data"," PC  Simulated Data","FA  Actual Data", " FA  Resampled Data", " FA  Simulated Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
       text.col = "green4", lty = c("solid","dotted","dashed","solid","dotted","dashed"),
       merge = TRUE, bg = 'gray90') 
       } else {legend("topright", c("PC  Actual Data", " PC  Resampled Data","FA  Actual Data", " FA  Resampled Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
       text.col = "green4", lty = c("solid","dotted","solid","dotted"),
       merge = TRUE, bg = 'gray90') 
   }
   
fa.test <- which(!(values > fam))[1]-1
pc.test <- which(!(e.values > pcm))[1] -1
cat("Parallel analysis suggests that ")
cat("the number of factors = ",fa.test, " and the number of components = ",pc.test,"\n")
 }
 #modified Oct 2, 2013 to plot (if desired) the simulated as well as the resampled data
 
 
 
 