\name{describe}
\alias{describe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Basic descriptive statistics useful for psychometrics }
\description{
There are many summary statistics available in R; this function
  provides the ones most useful for scale construction and item analysis in classic psychometrics. 
  Range is most useful for the first pass in a data set, to check for coding errors. 
}
\usage{
describe(x, digits = 2, na.rm = TRUE, skew = TRUE, ranges = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data frame or matrix}
  \item{digits}{ How many significant digits to report}
  \item{na.rm}{The default is to delete missing data }
  \item{skew}{ Should the skew and kurtosis be calculated? }
  \item{ranges}{ Should the range be calculated? }
}
\details{In basic data analysis it is vital to get basic descriptive statistics. 
Procedures such as \code{\link{summary}} and hmisc::describe do so.  The describe function in the \code{\link{psych}} package is meant to produce the most frequently requested stats in psychometric and psychology studies, and to produce them in an easy to read data.frame. The results from describe can be used in graphics functions (e.g., \code{\link{error.crosses}}).

The range statistics (min, max, range) are most useful for data checking to detect coding errors, and should be found in early analyses of the data.  

Although describe will work on data frames as well as matrices, it is important to realize that for data frames, descriptive statistics will be reported only for those variables where this makes sense (i.e., not for factors or for alphanumeric data). 

In a typical study, one might read the data in from the clipboard (\code{\link{read.clipboard}}), show the splom plot of the correlations (\code{\link{pairs.panels}}), and then describe the data. 


}
\value{
 A data.frame of the relevant statistics: \cr
   item name \cr
   item number \cr
   number of valid cases\cr
   mean\cr
   standard deviation\cr
   median\cr
   mad: median absolute deviation (from the median) \cr
   minimum\cr
   maximum\cr
   skew\cr
   kurtosis\cr
   standard error\cr
}

\note{Describe uses either the mean or colMeans functions depending upon whether the data are a data.frame or a matrix. The mean function supplies means for the columns of a data.frame, but the overall mean for a matrix.  Mean will throw a warning for non-numeric data, but colMeans stops with non-numeric data. Thus, the describe function uses either mean (for data frames) or colMeans (for matrices).  This is true for skew and kurtosi as well.}
\author{ 
 
 \url{http://personality-project.org/revelle.html} \cr

Maintainer: William Revelle \email{revelle@northwestern.edu} \cr
}
\seealso{ \code{\link{describe.by}}, \code{\link{skew}}, \code{\link{kurtosi}}, \code{\link{pairs.panels}}, \code{\link{read.clipboard}}, \code{\link{error.crosses}} }
\examples{

describe(attitude)
#           var  n  mean    sd median   mad min max range  skew kurtosis   se
#rating       1 30 64.63 12.17   65.5 10.38  40  85    45 -0.36    -0.77 2.22
#complaints   2 30 66.60 13.31   65.0 14.83  37  90    53 -0.22    -0.68 2.43
#privileges   3 30 53.13 12.24   51.5 10.38  30  83    53  0.38    -0.41 2.23
#learning     4 30 56.37 11.74   56.5 14.83  34  75    41 -0.05    -1.22 2.14
#raises       5 30 64.63 10.40   63.5 11.12  43  88    45  0.20    -0.60 1.90
#critical     6 30 74.77  9.89   77.5  7.41  49  92    43 -0.87     0.17 1.81
#advance      7 30 42.93 10.29   41.0  8.90  25  72    47  0.85     0.47 1.88
   
describe(attitude,skew=FALSE)   #attitude is taken from R data sets

#           var  n  mean    sd median   mad min max range   se
#rating       1 30 64.63 12.17   65.5 10.38  40  85    45 2.22
#complaints   2 30 66.60 13.31   65.0 14.83  37  90    53 2.43
#privileges   3 30 53.13 12.24   51.5 10.38  30  83    53 2.23
#learning     4 30 56.37 11.74   56.5 14.83  34  75    41 2.14
#raises       5 30 64.63 10.40   63.5 11.12  43  88    45 1.90
#critical     6 30 74.77  9.89   77.5  7.41  49  92    43 1.81
#advance      7 30 42.93 10.29   41.0  8.90  25  72    47 1.88





## The function is currently defined as

function (x, digits = 2,na.rm=TRUE,skew=TRUE,ranges=TRUE)   #basic stats after dropping non-numeric data
                                                 #slightly faster if we don't do skews
{                         #first, define a local function
    valid <- function(x) {sum(!is.na(x))}
   		 
    if (is.vector(x) )  {        #do it for vectors or 
    	    len  <- 1
    	    stats = matrix(rep(NA,8),ncol=8)    #create a temporary array
			stats[1, 1] <-  valid(x )			
			stats[1, 2] <-  mean(x, na.rm=na.rm )
			stats[1, 3] <-  median(x,na.rm=na.rm  )
			stats[1, 4] <-  min(x, na.rm=na.rm )
			stats[1, 5] <-  max(x, na.rm=na.rm )
			stats[1, 6] <-  skew(x,na.rm=na.rm  )
			stats[1,7] <-  mad(x,na.rm=na.rm) 
			stats[1,8] <-  kurtosi(x,na.rm=na.rm) 
			
    	}   else  {
    	len = dim(x)[2]     #do it for matrices or data.frames 
    	
   stats = matrix(rep(NA,len*8),ncol=8)    #create a temporary array
   stats[,1] <- apply(x,2,valid)
   if(is.matrix(x))  {stats[,2] <- colMeans(x, na.rm=na.rm )} else {stats[,2] <- mean(x,na.rm=na.rm)}
    if (skew) {stats[, 6] <-  skew(x,na.rm=na.rm  )
               stats[,8] <- kurtosi(x,na.rm=na.rm)}
    
    for (i in 1:len) {
    	if (is.numeric(x[,i])) {   #just do this for numeric data
			
		if (ranges) {
			stats[i, 3] <-  median(x[,i],na.rm=na.rm  )	
			stats[i,7] <-   mad(x[,i], na.rm=na.rm)
			stats[i, 4] <-  min(x[,i], na.rm=na.rm )
			stats[i, 5] <-  max(x[,i], na.rm=na.rm )
					} #ranges
		   		 }#is.numeric
        	}# i loop	
    	} #else loop
    if (ranges)
    	{if(skew){temp <-  data.frame(var = seq(1:len),n = stats[,1],mean=stats[,2], sd = sd(x,na.rm=TRUE), median = stats[, 
        3],mad = stats[,7], min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],skew = stats[, 6], kurtosis = stats[,8])}
         
      	 else {temp <-  data.frame(var = seq(1:len),n = stats[,1],mean=stats[,2], sd = sd(x,na.rm=TRUE), median = stats[, 
        3],mad = stats[,7],min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4])}}
        
        else {if(skew){temp <-  data.frame(var = seq(1:len),n = stats[,1],mean=stats[,2], sd = sd(x,na.rm=TRUE),skew = stats[, 6], kurtosis = stats[,8])}
       else {temp <-  data.frame(var = seq(1:len),n = stats[,1],mean=stats[,2], sd = sd(x,na.rm=TRUE))}}
                
    answer <-  round(data.frame(temp, se = temp$sd/sqrt(temp$n)),  digits)
     return(answer)
}



}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
