"fa.parallel" <-
function(x,n.obs=NULL,fa="both",main="Parallel Analysis Scree Plots",ntrials=20,error.bars=FALSE)  { 
	
	ci <- 1.96
	arrow.len <- .05
 nsub <- dim(x)[1]
 nvariables <- dim(x)[2]
 if (!is.null(n.obs)) { nsub <- n.obs 
  	rx <- x
  	if(dim(x)[1] != dim(x)[2]) {warning("You specified the number of subjects, implying a correlation matrix, but do not have a correlation matrix, correlations found ")
  	rx <- cor(x,use="pairwise") }   	 } else {

  	rx <- cor(x,use="pairwise")
 	}
 	
 	 if (nsub==nvariables) {warning("It seems as if you are using a correlation matrix, but have not specified the number of cases. The data are treated as if you had raw data.  ")}
  				
   valuesx  <- eigen(rx)$values 
   fa.valuesx  <- factor.pa(rx,warnings=FALSE)$values
 
  temp <- list(samp =list(),samp.fa <- list(),sim=list(),sim.fa=list())
   
 for (trials in 1:ntrials) {
   
    if(is.null(n.obs)) {sampledata <- matrix(sample(unlist(x),size=nsub*nvariables,replace=TRUE),nrow=nsub,ncol=nvariables) 
   					values.samp <- eigen(cor(sampledata,use="pairwise"))$values
   					temp[["samp"]][[trials]] <- values.samp
   					if (fa!= "pc") {
   						temp[["samp.fa"]][[trials]]  <- factor.pa(cor(sampledata,use="pairwise"),SMC="FALSE",warnings=FALSE)$values
          					}
                  } 
  
   simdata=matrix(rnorm(nsub*nvariables),nrow=nsub,ncol=nvariables)    #make up simulated data
   temp[["sim"]][[trials]] <- eigen(cor(simdata))$values
    
   if (fa!="pc") {
   
   		fa.values.sim <- factor.pa(cor(simdata),SMC="FALSE",warnings=FALSE)$values
   		 temp[["sim.fa"]][[trials]] <- 	fa.values.sim
}
   	}
   	
   if (fa!="pc") {ylabel <- "eigenvalues of principal components and factor analysis"} else { ylabel  <- "eigen values of principal components"}
   
   values.sim <- describe(t(matrix(unlist(temp[["sim"]]),ncol=ntrials)))
  if(is.null(n.obs)) { values.samp <- describe(t(matrix(unlist(temp[["samp"]]),ncol=ntrials)))}
  
    ymax <- max(valuesx,values.sim$mean)
   
if (fa !="fa") {plot(valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=4) 
	points(values.sim$mean,type ="l",lty="dotted",pch=4)
	if(error.bars) {
      for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- values.sim$mean[i]
         yse <-  values.sim$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	
	
	if(is.null(n.obs)) {points(values.samp$mean,type ="l",lty="dashed",pch=4)
	if(error.bars) {
      for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- values.samp$mean[i]
         yse <-  values.samp$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	       
	           }
	}

if (fa !="pc" ) { if (fa=="fa") { ylabel <-  "eigen values of principal factors"
            plot(fa.valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=4) }
            fa.values.sim <- describe(t(matrix(unlist(temp[["sim.fa"]]),ncol=ntrials)))
        	points(fa.values.sim$mean,type ="l",lty="dotted",pch=2)
        if(error.bars) {
         for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- fa.values.sim$mean[i]
         yse <-  fa.values.sim$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	
        	  
	if(is.null(n.obs)) {fa.values.samp <- describe(t(matrix(unlist(temp[["samp.fa"]]),ncol=ntrials)))
						points(fa.values.samp$mean,type ="l",lty="dashed",pch=2)
						 if(error.bars) {
         for (i in 1:dim(values.sim)[1])  
    	{
    	 ycen <- fa.values.samp$mean[i]
         yse <-  fa.values.samp$se[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} }
    	} 
					if (fa !="fa") 	points(fa.valuesx,type ="b",lty="solid",pch=2)
	points(fa.values.sim,type ="l",lty="dotted",pch=2)
	if(is.null(n.obs)) {points(fa.values.samp$mean,type ="l",lty="dashed",pch=2)}
        } 


if(is.null(n.obs)) {
legend("topright", c("  PC  Actual Data", "  PC  Simulated Data", " PC  Resampled Data","  FA  Actual Data", "  FA  Simulated Data", " FA  Resampled Data"), col = c(3,4,5,6,7,8),pch=c(4,NA,NA,2,NA,NA),
       text.col = "green4", lty = c("solid","dotted", "dashed","solid","dotted", "dashed"),
       merge = TRUE, bg = 'gray90')} else {

       legend("topright", c("PC  Actual Data", " PC  Simulated Data","FA  Actual Data", " FA  Simulated Data"), col = c(3,4,5,6),pch=c(4,NA,2,NA),
       text.col = "green4", lty = c("solid","dotted","solid","dotted"),
       merge = TRUE, bg = 'gray90')}
       
abline(h=1)
if (fa!="pc") {abline(h=0) }}

 