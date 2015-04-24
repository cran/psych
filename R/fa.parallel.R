#1/2/14  switched the n.iter loop to a mclapply loop to allow for multicore parallel processing

"fa.parallel" <-
function(x,n.obs=NULL,fm="minres",fa="both",main="Parallel Analysis Scree Plots",n.iter=20,error.bars=FALSE,se.bars=TRUE,SMC=FALSE,ylabel=NULL,show.legend=TRUE,sim=TRUE,quant=.95,cor="cor",use="pairwise")  { 
 cl <- match.call()
# if(!require(parallel)) {message("The parallel package needs to be installed to run mclapply")}
	
	ci <- 1.96
	arrow.len <- .05
 nsub <- dim(x)[1]
 nvariables <- dim(x)[2]
 resample <- TRUE  #this is used as a flag for correlation matrices
 if((nsub == nvariables) && !sim)  {warning("You specified a correlation matrix, but asked to just resample (sim was set to FALSE).  This is impossible, so sim is set to TRUE")
  	  sim <- TRUE
  	  resample <- FALSE}
 if (!is.null(n.obs)) { nsub <- n.obs 
  	rx <- x
  	resample <- FALSE
  	if(dim(x)[1] != dim(x)[2]) {warning("You specified the number of subjects, implying a correlation matrix, but do not have a correlation matrix, correlations found ")
  #	rx <- cor(x,use="pairwise") 
  #add the option to choose the type of correlation, this allows us to do fa.parallel.poly inside fa.parallel
  switch(cor, 
       cor = {rx <- cor(x,use=use)},
       cov = {rx <- cov(x,use=use) 
              covar <- TRUE},
       tet = {rx <- tetrachoric(x)$rho},
       poly = {rx <- polychoric(x)$rho},
       mixed = {rx <- mixed.cor(x,use=use)$rho},
       Yuleb = {rx <- YuleCor(x,,bonett=TRUE)$rho},
       YuleQ = {rx <- YuleCor(x,1)$rho},
       YuleY = {rx <- YuleCor(x,.5)$rho } 
       )
  	
  	if(!sim) {warning("You specified a correlation matrix, but asked to just resample (sim was set to FALSE).  This is impossible, so sim is set to TRUE")
  	  sim <- TRUE
  	  resample <- FALSE}
  	}   	 } else {
  	if (nsub==nvariables) {warning("It seems as if you are using a correlation matrix, but have not specified the number of cases. The number of subjects is arbitrarily set to be 100  ") 
  	rx <- x
  	nsub = 100
  	n.obs=100
  	resample <- FALSE
      }  else {
  
  	switch(cor, 
       cor = {rx <- cor(x,use=use)},
       cov = {rx <- cov(x,use=use) 
              covar <- TRUE},
       tet = {rx <- tetrachoric(x)$rho},
       poly = {rx <- polychoric(x)$rho},
       mixed = {rx <- mixed.cor(x,use=use)$rho},
       Yuleb = {rx <- YuleCor(x,,bonett=TRUE)$rho},
       YuleQ = {rx <- YuleCor(x,1)$rho},
       YuleY = {rx <- YuleCor(x,.5)$rho } 
       )
 	} }
 	
 	 
  				
   valuesx  <- eigen(rx)$values #these are the PC values
   if(SMC) {diag(rx) <- smc(rx)
   fa.valuesx <- eigen(rx)$values} else {
   fa.valuesx  <- fa(rx,fm=fm,warnings=FALSE)$values}  #these are the FA values
 
  temp <- list(samp =vector("list",n.iter),samp.fa = vector("list",n.iter),sim=vector("list",n.iter),sim.fa=vector("list",n.iter))
   
#parallel processing starts here  - the more cores the better!
   templist <- mclapply(1:n.iter,function(XX) {
    if(is.null(n.obs)) {
    #Sample the data, column wise (to keep the basic distributional properties, but making the correlations 0 (on average)
    bad <- TRUE
    while(bad) {sampledata <- matrix(apply(x,2,function(y) sample(y,nsub,replace=TRUE)),ncol=nvariables) #do it column wise
                    
    switch(cor,          #we can do a number of different types of correlations
       cor = {C <- cor(sampledata,use=use)},
       cov = {C <- cov(sampledata,use=use) 
              covar <- TRUE},
       tet = {C <- tetrachoric(sampledata)$rho},
       poly = {C <- polychoric(sampledata)$rho},
       mixed = {C <- mixed.cor(sampledata,use=use)$rho},
       Yuleb = {C <- YuleCor(sampledata,,bonett=TRUE)$rho},
       YuleQ = {C <- YuleCor(sampledata,1)$rho},
       YuleY = {C <- YuleCor(sampledata,.5)$rho } 
       )  
                    bad <- any(is.na(C))   #some (not frequently) correlations will be improper, particularly if sampling from sparse matrices                 
           }  #Try resampling until we get a correlation matrix that works                    
   					values.samp <- eigen(C)$values
   					temp[["samp"]] <- values.samp
   					if (fa!= "pc") {
   					 if(SMC) {sampler <- C 
   					          diag(sampler) <- smc(sampler)
   					 temp[["samp.fa"]]<- eigen(sampler)$values} else {
   					 temp[["samp.fa"]]  <- fa(C,fm=fm,SMC=FALSE,warnings=FALSE)$values
          					}
          					}
                  } 
  
  if(sim) { simdata=matrix(rnorm(nsub*nvariables),nrow=nsub,ncol=nvariables)    #create simulated data based upon normal theory
   sim.cor <- cor(simdata)   #we must use correlations based upon Pearson here, because we are simulating the data
  
   temp[["sim"]] <- eigen(sim.cor)$values
    
   if (fa!="pc") {
        if(SMC) { diag(sim.cor) <- smc(sim.cor)
   					 temp[["sim.fa"]]<- eigen(sim.cor)$values} else {fa.values.sim <- fa(sim.cor,fm=fm,SMC=FALSE,warnings=FALSE)$values
   		 temp[["sim.fa"]]    <- fa.values.sim
}}}
   replicates <- list(samp=temp[["samp"]],samp.fa=temp[["samp.fa"]],sim=temp[["sim"]],sim.fa=temp[["sim.fa"]])
   	})
#parallelism stops here
#now combine the results   	
   	
   if(is.null(ylabel)) {if (fa!="pc") {ylabel <- "eigenvalues of principal components and factor analysis"} else { ylabel  <- "eigen values of principal components"}}
   values<- t(matrix(unlist(templist),ncol=n.iter))
   
  
   values.sim.mean=colMeans(values,na.rm=TRUE)
    if(!missing(quant)) {values.ci = apply(values,2,function(x) quantile(x,quant))} else {values.ci <- values.sim.mean}
   if(se.bars) {values.sim.se <- apply(values,2,sd,na.rm=TRUE)/sqrt(n.iter)} else {values.sim.se <- apply(values,2,sd,na.rm=TRUE)}
  
    ymax <- max(valuesx,values.sim.mean)
    
switch(fa, 
  pc = { plot(valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Component Number",pch=4,col="blue")
        if(resample) { sim.pcr <- values.sim.mean[1:nvariables] 
                       sim.pcr.ci <- values.ci[1:nvariables]
                       sim.se.pcr <- values.sim.se[1:nvariables]
         points(sim.pcr,type ="l",lty="dashed",pch=4,col="red")} else {sim.pcr <- NA
              sim.se.pc <- NA}
       if(sim) {    
        if(resample) {sim.pc <- values.sim.mean[(nvariables+1):(2*nvariables)] 
                      sim.pc.ci <- values.ci[(nvariables+1):(2*nvariables)] 
                      sim.se.pc <- values.sim.se[(nvariables+1):(2*nvariables)] 
        } else {sim.pc <- values.sim.mean[1:nvariables]
                sim.pc.ci <- values.ci[1:nvariables]
                sim.se.pc <- values.sim.se[1:nvariables]}
        points(sim.pc,type ="l",lty="dotted",pch=4,col="red")
        pc.test <- which(!(valuesx > sim.pc.ci))[1]-1} else { 
           sim.pc <- NA
           sim.pc.ci <- NA 
           sim.se.pc <- NA
           pc.test <- which(!(valuesx > sim.pcr.ci))[1]-1  }
        fa.test <- NA
        sim.far <- NA
        sim.fa <- NA 
    },
  fa = { ylabel <-  "eigen values of principal factors"
        plot(fa.valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=2,col="blue")
         sim.se.pc <- NA
          if(resample) {sim.far <- values.sim.mean[(nvariables+1):(2*nvariables)]
                        sim.far.ci <- values.ci[(nvariables+1):(2*nvariables)]
                        sim.se.far <- values.sim.se[(nvariables+1):(2*nvariables)]
             points(sim.far,type ="l",lty="dashed",pch=2,col="red")}
        if(sim) { if(resample) {sim.fa <- values.sim.mean[(3*nvariables+1):(4*nvariables)] 
                                sim.fa.ci <- values.ci[(3*nvariables+1):(4*nvariables)]
                                sim.se.fa <- values.sim.se[(3*nvariables+1):(4*nvariables)] } else {
                 sim.fa <- values.sim.mean[(nvariables+1):(2*nvariables)]
                 sim.fa.ci <- values.sim.mean[(nvariables+1):(2*nvariables)]
                sim.se.fa <- values.sim.se[(nvariables+1):(2*nvariables)]}
        	points(sim.fa,type ="l",lty="dotted",pch=2,col="red")
        	fa.test <- which(!(fa.valuesx > sim.fa.ci))[1]-1 
        } else {sim.fa <- NA 
        fa.test <- which(!(fa.valuesx > sim.far.ci))[1]-1 }
         sim.pc <- NA
         sim.pcr <- NA
         sim.se.pc <- NA
         pc.test <- NA },
  both = { plot(valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor/Component Number",pch=4,col="blue")
           points(fa.valuesx,type="b",pch=2,col="blue")
     if(sim) {  
   	 if(resample) {
   	 sim.pcr <- values.sim.mean[1:nvariables]
   	 sim.pcr.ci <- values.ci[1:nvariables] 
   	 sim.se.pcr <- values.sim.se[1:nvariables]  
   	 sim.far <- values.sim.mean[(nvariables+1):(2*nvariables)]
   	 sim.se.far <- values.sim.se[(nvariables+1):(2*nvariables)]
   	 sim.far.ci <- values.ci[(nvariables+1):(2*nvariables)]
     sim.pc <- values.sim.mean[(2*nvariables+1):(3*nvariables)]
     sim.pc.ci <- values.ci[(2*nvariables+1):(3*nvariables)]
     sim.se.pc <- values.sim.se[(2*nvariables+1):(3*nvariables)]
   	 sim.fa <- values.sim.mean[(3*nvariables+1):(4*nvariables)]
   	 sim.fa.ci <- values.ci[(3*nvariables+1):(4*nvariables)]
   	 sim.se.fa <- values.sim.se[(3*nvariables+1):(4*nvariables)]
   	 pc.test <- which(!(valuesx > sim.pcr.ci))[1]-1  
     fa.test <- which(!(fa.valuesx > sim.far.ci))[1]-1  } else { #not resampling, just sim
   	 sim.pc <- values.sim.mean[1:nvariables] 
   	 sim.pc.ci <- values.ci[1:nvariables] 
   	 sim.se.pc <- values.sim.se[1:nvariables]  
   	 sim.fa <- values.sim.mean[(nvariables+1):(2*nvariables)]
   	 sim.fa.ci <- values.ci[(nvariables+1):(2*nvariables)]
   	 sim.se.fa <- values.sim.se[(nvariables+1):(2*nvariables)]
   	  pc.test <- which(!(valuesx > sim.pc.ci))[1]-1  
      fa.test <- which(!(fa.valuesx > sim.fa.ci))[1]-1 }
   	 
     points(sim.pc,type ="l",lty="dotted",pch=4,col="red")
     points(sim.fa,type ="l",lty="dotted",pch=4,col="red")
     sim.pcr <- sim.far <- NA
#points(sim.pcr,type ="l",lty="dashed",pch=4,col="red")
#     points(sim.far,type ="l",lty="dashed",pch=4,col="red")
      pc.test <- which(!(valuesx > sim.pc.ci))[1]-1 
      fa.test <- which(!(fa.valuesx > sim.fa.ci))[1]-1  
     } else {      #sim is false
    sim.pcr <- values.sim.mean[1:nvariables]
    sim.pcr.ci <- values.ci[1:nvariables]
    sim.se.pcr <- values.sim.se[1:nvariables]
    sim.far <- values.sim.mean[(nvariables+1):(2*nvariables)]
     sim.far.ci <- values.ci[(nvariables+1):(2*nvariables)]
    sim.se.far <- values.sim.se[(nvariables+1):(2*nvariables)]
    sim.fa <- NA
    sim.pc <- NA
    sim.se.fa <- NA
    sim.se.pc <- NA
    pc.test <- which(!(valuesx > sim.pcr.ci))[1]-1  
    fa.test <- which(!(fa.valuesx > sim.far.ci))[1]-1 } 
    if(resample) { points(sim.pcr,type ="l",lty="dashed",pch=4,col="red")
        points(sim.far,type ="l",lty="dashed",pch=4,col="red") 
     } } 
)  

 

 	if(error.bars) {  
 	if(!any(is.na(sim.pc))) {  
       for (i in 1:length(sim.pc))  {
    	 ycen <- sim.pc[i]
         yse <-  sim.se.pc[i]
    	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} 
    	 }
  if(!any(is.na(sim.pcr))) {    
       for (i in 1:length(sim.pcr))   {
     	 ycen <- sim.pcr[i]
          yse <-  sim.se.pcr[i]
     	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} 
    	}
  if(!any(is.na(sim.fa))) {      	       
 	   for (i in 1:length(sim.fa))  {
     	 ycen <- sim.fa[i]
          yse <-  sim.se.fa[i]
     	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} 
     	}
 if(!any(is.na(sim.far))) {      	
          for (i in 1:length(sim.far))  {
     	 ycen <- sim.far[i]
          yse <-  sim.se.far[i]
     	 arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} 
    	} 
    	}


if(show.legend) {
if(is.null(n.obs)) {   #that is, do we have real data or a correlation matrix
switch(fa,  
both = {if(sim) {legend("topright", c("  PC  Actual Data", "  PC  Simulated Data", " PC  Resampled Data","  FA  Actual Data", "  FA  Simulated Data", " FA  Resampled Data"),
       col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
       text.col = "green4", lty = c("solid","dotted", "dashed","solid","dotted", "dashed"),
      merge = TRUE, bg = 'gray90')
      } else {legend("topright", c("  PC  Actual Data",  " PC  Resampled Data","  FA  Actual Data",  " FA  Resampled Data"),
       col = c("blue","red","blue","red"),pch=c(4,NA,2,NA,NA),
       text.col = "green4", lty = c("solid","dashed", "solid","dashed"),
      merge = TRUE, bg = 'gray90')}},
       
pc = {if(sim) {legend("topright", c("  PC  Actual Data", "  PC  Simulated Data", " PC  Resampled Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
       text.col = "green4", lty = c("solid","dotted", "dashed","solid","dotted", "dashed"),
       merge = TRUE, bg = 'gray90')} else {
       legend("topright", c("  PC  Actual Data",  " PC  Resampled Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
       text.col = "green4", lty = c("solid", "dashed","solid","dotted", "dashed"),
       merge = TRUE, bg = 'gray90')
       } } ,  
       
fa = {if(sim) {legend("topright", c("  FA  Actual Data", "  FA  Simulated Data", " FA  Resampled Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
       text.col = "green4", lty = c("solid","dotted", "dashed","solid","dotted", "dashed"),
       merge = TRUE, bg = 'gray90')} else {
        legend("topright", c("  FA  Actual Data",  " FA  Resampled Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
       text.col = "green4", lty = c("solid", "dashed","solid","dotted", "dashed"),
       merge = TRUE, bg = 'gray90')
       }
       }   
       ) } else {
switch(fa,

 both= {      legend("topright", c("PC  Actual Data", " PC  Simulated Data","FA  Actual Data", " FA  Simulated Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
       text.col = "green4", lty = c("solid","dotted","solid","dotted"),
       merge = TRUE, bg = 'gray90')},
       
  pc= {   legend("topright", c("PC  Actual Data", " PC  Simulated Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
       text.col = "green4", lty = c("solid","dotted","solid","dotted"),
       merge = TRUE, bg = 'gray90')},
fa =  {legend("topright", c("FA  Actual Data", " FA  Simulated Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
       text.col = "green4", lty = c("solid","dotted","solid","dotted"),
       merge = TRUE, bg = 'gray90')})}
   }
   

colnames(values) <- paste0("Sim",1:ncol(values))
if(fa!= "pc") abline(h=1)
results <- list(fa.values = fa.valuesx,pc.values=valuesx,pc.sim=sim.pc,pc.simr = sim.pcr,fa.sim=sim.fa,fa.simr = sim.far,nfact=fa.test,ncomp=pc.test, Call=cl) 



if (fa == "pc" )  {
      colnames(values)[1:nvariables] <- paste0("C",1:nvariables)
    } else {

 
colnames(values)[1:(2*nvariables)] <- c(paste0("C",1:nvariables),paste0("F",1:nvariables))
if(sim) {
if(resample) colnames(values)[(2*nvariables +1):ncol(values)] <- c(paste0("CSim",1:nvariables),paste0("Fsim",1:nvariables))
} 
results$nfact <- fa.test}

results$ncomp <- pc.test
results$values <- values
cat("Parallel analysis suggests that ")
cat("the number of factors = ",fa.test, " and the number of components = ",pc.test,"\n")
class(results) <- c("psych","parallel")
return(invisible(results))
}


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
 #modified June 09, 2013 to fix the case for just PC tests
 #modified October 2, 2013 to sample each column of data separately.  
 #modified March 23, 2014 to check for bad resamples in case of very sparse data
 #also modified to just resample if desired
 
"test.fa.parallel" <- function(x,n.obs=NULL) {
fp <- list()

fp[[1]] <- fa.parallel(x,n.obs=n.obs)
fp[[2]] <-fa.parallel(x,sim=FALSE,n.obs=n.obs)
fp[[3]] <- fa.parallel(x,error.bars=TRUE,n.obs=n.obs)
fp[[4]] <- fa.parallel(x,error.bars=TRUE,sim=FALSE,n.obs=n.obs)
fp[[5]] <- fa.parallel(x,fa="fa",n.obs=n.obs)
fp[[6]] <- fa.parallel(x,fa="fa",sim=FALSE,n.obs=n.obs)
fp[[7]] <- fa.parallel(x,fa="fa",error.bars=TRUE,n.obs=n.obs)
fp[[8]] <- fa.parallel(x,fa="fa",error.bars=TRUE,sim=FALSE,n.obs=n.obs)
fp[[9]] <- fa.parallel(x,fa="pc",n.obs=n.obs)
fp[[10]] <- fa.parallel(x,fa="pc",sim=FALSE,n.obs=n.obs)
fp[[11]] <- fa.parallel(x,fa="pc",error.bars=TRUE,n.obs=n.obs)
fp[[12]] <- fa.parallel(x,fa="pc",error.bars=TRUE,sim=FALSE,n.obs=n.obs)
fp[[13]] <- fa.parallel(x,quant=.95,n.obs=n.obs)
return(fp)
}

print.psych.parallel <- function(x,digits=2) {
cat("Call: ")
              print(x$Call) 
              if(!is.null(x$fa.values) & !is.null(x$pc.values) ) {
                  parallel.df <- data.frame(fa=x$fa.values,fa.sim =x$fa.sim,pc= x$pc.values,pc.sim =x$pc.sim)
		fa.test <- x$nfact
		pc.test <- x$ncomp
		
		cat("Parallel analysis suggests that ")
		cat("the number of factors = ",fa.test, " and the number of components = ",pc.test,"\n")
                  cat("\n Eigen Values of \n")
                  
                  colnames(parallel.df) <- c("Original factors","Simulated data","Original components", "simulated data")}
              if(is.na(fa.test) ) fa.test <- 0
              if(is.na(pc.test)) pc.test <- 0
               if(!any(is.na(parallel.df)))  {print(round(parallel.df[1:max(fa.test,pc.test),],digits))}  else {
              if(!is.null(x$fa.values)) {cat("\n eigen values of factors\n")
              print(round(x$fa.values,digits))}
               if(!is.null(x$fa.sim)){cat("\n eigen values of simulated factors\n") 
               print(round(x$fa.sim,digits))}
              if(!is.null(x$pc.values)){cat("\n eigen values of components \n")
              print(round(x$pc.values,digits))}
             
              if(!is.null(x$pc.sim)) {cat("\n eigen values of simulated components\n") 
              print(round(x$pc.sim,digits=digits))}
            }  

}




 