"fa.parallel" <-
function(x,n.obs=NULL,fa="both",main="Parallel Analysis Scree Plots")  {     
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
   						 fa.values.samp <- factor.pa(cor(sampledata,use="pairwise"))$values
          					}
                  } 
    
   simdata=matrix(rnorm(nsub*nvariables),nrow=nsub,ncol=nvariables)    #make up simulated data
   values.sim <- eigen(cor(simdata))$values

   valuesx  <- eigen(rx)$values 
   ylabel  <- "eigen values of principal components"
   ymax <- max(valuesx)
   
   if (fa!="pc") {
   
   	fa.values.sim <- factor.pa(cor(simdata))$values
   	fa.valuesx  <- factor.pa(rx)$values
    ylabel <-  "eigen values of principal components and factors"
   	}
    

if (fa !="fa") {plot(valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=4) 
	points(values.sim,type ="b",lty="dotted",pch=4)
	if(is.null(n.obs)) {points(values.samp,type ="b",lty="dashed",pch=4)} }


if (fa !="pc" ) { if (fa=="fa") {
             ylabel <-  "eigen values of principal factors"
             plot(fa.valuesx,type="b", main = main,ylab=ylabel ,ylim=c(0,ymax),xlab="Factor Number",pch=4) 
        	points(fa.values.sim,type ="b",lty="dotted",pch=2)
	if(is.null(n.obs)) {points(fa.values.samp,type ="b",lty="dashed",pch=2)}} else
	points(fa.valuesx,type ="b",lty="solid",pch=2)
	points(fa.values.sim,type ="b",lty="dotted",pch=2)
	if(is.null(n.obs)) {points(fa.values.samp,type ="b",lty="dashed",pch=2)}
        } 


if(is.null(n.obs)) {
legend("topright", c("  PC  Actual Data", "  PC  Simulated Data", " PC  Resampled Data","  FA  Actual Data", "  FA  Simulated Data", " FA  Resampled Data"), col = c(3,4,5,6,7,8),pch=c(4,4,4,2,2,2),
       text.col = "green4", lty = c("solid","dotted", "dashed","solid","dotted", "dashed"),
       merge = TRUE, bg = 'gray90')} else {

       legend("topright", c("PC  Actual Data", " PC  Simulated Data","FA  Actual Data", " FA  Simulated Data"), col = c(3,4,5,6),pch=c(4,4,2,2),
       text.col = "green4", lty = c("solid","dotted","solid","dotted"),
       merge = TRUE, bg = 'gray90')}
       
abline(h=1)
if (fa!="pc") {abline(h=0) }}

 