"fa.parallel" <-
function(x,n.obs=NULL,fm="minres",fa="both",main="Parallel Analysis Scree Plots",n.iter=20,error.bars=FALSE,SMC=FALSE,ylabel=NULL,show.legend=TRUE)  { 
 cl <- match.call()
	
	ci <- 1.96
	arrow.len <- .05
 nsub <- dim(x)[1]
 nvariables <- dim(x)[2]
 if (!is.null(n.obs)) { nsub <- n.obs 
  	rx <- x
  	if(dim(x)[1] != dim(x)[2]) {warning("You specified the number of subjects, implying a correlation matrix, but do not have a correlation matrix, correlations found ")
  	rx <- cor(x,use="pairwise") }   	 } else {
  	if (nsub==nvariables) {warning("It seems as if you are using a correlation matrix, but have not specified the number of cases. The number of subjects is arbitrarily set to be 100  ") 
  	rx <- x
  	nsub = 100
  	n.obs=100}  else {
  	rx <- cor(x,use="pairwise")
 	} }
 	
 	 
  				
   valuesx  <- eigen(rx)$values #these are the PC values
   if(SMC) {diag(rx) <- smc(rx)
   fa.valuesx <- eigen(rx)$values} else {
   fa.valuesx  <- fa(rx,fm=fm,warnings=FALSE)$values}  #these are the FA values
 
  temp <- list(samp =list(),samp.fa <- list(),sim=list(),sim.fa=list())
   
 for (trials in 1:n.iter) {
   
    if(is.null(n.obs)) {sampledata <- matrix(sample(unlist(x),size=nsub*nvariables,replace=TRUE),nrow=nsub,ncol=nvariables) 
   					values.samp <- eigen(cor(sampledata,use="pairwise"))$values
   					temp[["samp"]][[trials]] <- values.samp
   					if (fa!= "pc") {
   					 if(SMC) {sampler <- cor(sampledata,use="pairwise")
   					          diag(sampler) <- smc(sampler)
   					 temp[["samp.fa"]][[trials]] <- eigen(sampler)$values} else {
   						temp[["samp.fa"]][[trials]]  <- fa(cor(sampledata,use="pairwise"),fm=fm,SMC=FALSE,warnings=FALSE)$values
          					}
          					}
                  } 
  
   simdata=matrix(rnorm(nsub*nvariables),nrow=nsub,ncol=nvariables)    #make up simulated data
   sim.cor <- cor(simdata)
   temp[["sim"]][[trials]] <- eigen(sim.cor)$values
    
   if (fa!="pc") {
        if(SMC) { diag(sim.cor) <- smc(sim.cor)
   					 temp[["sim.fa"]][[trials]] <- eigen(sim.cor)$values} else 
   		{fa.values.sim <- fa(sim.cor,fm=fm,SMC=FALSE,warnings=FALSE)$values
   		 temp[["sim.fa"]][[trials]] <- 	fa.values.sim
}}
   	}
   	
   	
   	
   if(is.null(ylabel)) {if (fa!="pc") {ylabel <- "eigenvalues of principal components and factor analysis"} else { ylabel  <- "eigen values of principal components"}}
   
   values.sim <- describe(t(matrix(unlist(temp[["sim"]]),ncol=n.iter)))
  if(is.null(n.obs)) { values.samp <- describe(t(matrix(unlist(temp[["samp"]]),ncol=n.iter)))}
  
    ymax <- max(valuesx,values.sim$mean)
   
if (fa !="fa") {plot(valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=4,col="blue") 
	points(values.sim$mean,type ="l",lty="dotted",pch=4,col="red")
	if(error.bars) {
      for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- values.sim$mean[i]
         yse <-  values.sim$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	
	
	if(is.null(n.obs)) {points(values.samp$mean,type ="l",lty="dashed",pch=4,col="red")
	if(error.bars) {
      for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- values.samp$mean[i]
         yse <-  values.samp$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	       
	           }
	}

if (fa !="pc" ) { if (fa=="fa") { ylabel <-  "eigen values of principal factors"
            plot(fa.valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=4,col="blue") }
            fa.values.sim <- describe(t(matrix(unlist(temp[["sim.fa"]]),ncol=n.iter)))
        	points(fa.values.sim$mean,type ="l",lty="dotted",pch=2,col="red")
        if(error.bars) {
         for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- fa.values.sim$mean[i]
         yse <-  fa.values.sim$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	
        	  
	if(is.null(n.obs)) {fa.values.samp <- describe(t(matrix(unlist(temp[["samp.fa"]]),ncol=n.iter)))
						points(fa.values.samp$mean,type ="l",lty="dashed",pch=2,col="red")
						 if(error.bars) {
         for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- fa.values.samp$mean[i]
         yse <-  fa.values.samp$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	} 
					if (fa !="fa") 	points(fa.valuesx,type ="b",lty="solid",pch=2,col="blue")
	points(fa.values.sim,type ="l",lty="dotted",pch=2,col="red")
	if(is.null(n.obs)) {points(fa.values.samp$mean,type ="l",lty="dashed",pch=2,col="red")}
        } 

if(show.legend) {
if(is.null(n.obs)) {
legend("topright", c("  PC  Actual Data", "  PC  Simulated Data", " PC  Resampled Data","  FA  Actual Data", "  FA  Simulated Data", " FA  Resampled Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
       text.col = "green4", lty = c("solid","dotted", "dashed","solid","dotted", "dashed"),
       merge = TRUE, bg = 'gray90')} else {

       legend("topright", c("PC  Actual Data", " PC  Simulated Data","FA  Actual Data", " FA  Simulated Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
       text.col = "green4", lty = c("solid","dotted","solid","dotted"),
       merge = TRUE, bg = 'gray90')}
   }
abline(h=1)
if (fa!="pc") {abline(h=0) }
results <- list(fa.values = fa.valuesx,fa.sim = fa.values.sim,pc.values=valuesx,pc.sim=values.sim,Call=cl)
fa.test <- which(!(fa.valuesx > fa.values.sim$mean))[1]-1
pc.test <- which(!(valuesx > values.sim$mean))[1] -1
results$nfact <- fa.test
results$ncomp <- pc.test
cat("Parallel analysis suggests that ")
cat("the number of factors = ",fa.test, " and the number of components = ",pc.test,"\n")
class(results) <- c("psych","parallel")
return(invisible(results))}


#a cut down plotting function
"plot.fa.parallel" <- 
function(x,n.obs,fa,show.legend,error.bars,main="Parallel Analysis Scree Plots",...) {
if(missing(n.obs)) n.obs <- NULL
if(missing(fa)) fa <- "both"
if(missing(show.legend)) show.legend <- TRUE
if(missing(error.bars)) error.bars <- FALSE
ci <- 1.96
arrow.len <- .05
fa.valuesx <- x$fa.values
fa.values.sim <- x$fa.sim
valuesx <- x$pc.values
values.sim <- x$pc.sim
ymax <- max(valuesx,values.sim$mean)
ylabel <-  "eigen values of principal factors"
if (!is.null(valuesx)) {
            plot(valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=4,col="blue") }
           
        	points(values.sim$mean,type ="l",lty="dotted",pch=2,col="red")
        if(error.bars) {
         for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- fa.values.sim$mean[i]
         yse <-  fa.values.sim$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	
        	  
	
		points(fa.values.sim$mean,type ="l",lty="dashed",pch=2,col="red")
						 if(error.bars) {
         for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- fa.values.sim$mean[i]
         yse <-  fa.values.sim$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	
    	
					if (fa !="fa") 	points(fa.valuesx,type ="b",lty="solid",pch=2,col="blue")
	points(fa.values.sim,type ="l",lty="dotted",pch=2,col="red")
	if(is.null(n.obs)) {points(fa.values.sim$mean,type ="l",lty="dashed",pch=2,col="red")}
           

if(show.legend) {

       legend("topright", c("PC  Actual Data", " PC  Simulated Data","FA  Actual Data", " FA  Simulated Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
       text.col = "green4", lty = c("solid","dotted","solid","dotted"),
       merge = TRUE, bg = 'gray90')
   }
abline(h=1)
if (fa!="pc") {abline(h=0) }
}




 