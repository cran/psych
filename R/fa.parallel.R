"fa.parallel" <-
function(x,n.obs=1000,fa="both",main="Parallel Analysis Scree Plots",ntrials=20)  {     
 nsub <- dim(x)[1]
 nvariables <- dim(x)[2]
 if (!is.null(n.obs)) { nsub <- n.obs 
  	rx <- x
  	if(dim(x)[1] != dim(x)[2]) {warning("You specified the number of subjects, implying a correlation matrix, but do not have a correlation matrix, correlations found ")
  	rx <- cor(x,use="pairwise") }   	 } else {

  	rx <- cor(x,use="pairwise")
 	}
  
  if(is.null(n.obs)) {sampledata <- matrix(sample(unlist(x),size=nsub*nvariables,replace=TRUE),nrow=nsub,ncol=nvariables) 
   					values.samp <- eigen(cor(sampledata,use="pairwise"))$values
   					 if (nsub==nvariables) {warning("It seems as if you are using a correlation matrix, but have not specified the number of cases. The data are treated as if you had raw data.  ")}
  					if (fa!= "pc") {
   						 fa.values.samp <- factor.pa(cor(sampledata,use="pairwise"),warnings=FALSE)$values
          					}
                  } 
    
   temp <- rep(0,nvariables)
   temp.fa <- rep(0,nvariables)
   
   for (trials in 1:ntrials) {
   simdata=matrix(rnorm(nsub*nvariables),nrow=nsub,ncol=nvariables)    #make up simulated data
   values.sim <- eigen(cor(simdata))$values

   valuesx  <- eigen(rx)$values 
   temp <- temp + valuesx
  
  
   
   if (fa!="pc") {
   
   		fa.values.sim <- factor.pa(cor(simdata),SMC="FALSE",warnings=FALSE)$values
   	
   		fa.valuesx  <- factor.pa(rx,warnings=FALSE)$values
   		temp.fa <- temp.fa + fa.valuesx
   		ylabel <-  "eigen values of principal components and factors"
   				}
   	}
    ylabel  <- "eigen values of principal components"
   valuesx <- temp/ntrials
   fa.values.x <- temp.fa/ntrials
   ymax <- max(valuesx)
   

if (fa !="fa") {plot(valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=4) 
	points(values.sim,type ="l",lty="dotted",pch=4)
	if(is.null(n.obs)) {points(values.samp,type ="l",lty="dashed",pch=4)} }


if (fa !="pc" ) { if (fa=="fa") {
             ylabel <-  "eigen values of principal factors"
             plot(fa.valuesx,type="l", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=4) 
        	points(fa.values.sim,type ="l",lty="dotted",pch=2)
	if(is.null(n.obs)) {points(fa.values.samp,type ="l",lty="dashed",pch=2)}} else
	points(fa.valuesx,type ="b",lty="solid",pch=2)
	points(fa.values.sim,type ="l",lty="dotted",pch=2)
	if(is.null(n.obs)) {points(fa.values.samp,type ="l",lty="dashed",pch=2)}
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

 