\name{describe}
\alias{describe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Basic descriptive statistics useful for psychometrics }
\description{
There are many summary statistics available in R; this function
  provides the ones most useful for scale construction and item analysis in classic psychometrics.
  Although describe (optionally) will calculate skew, for large data sets this produces a noticable slowing. 
  Range is most useful for the first pass in a data set, to check for coding errors. 
}
\usage{
describe(x, digits = 2, na.rm = TRUE, skew = FALSE, ranges = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data frame or matrix}
  \item{digits}{ How many significant digits to report}
  \item{na.rm}{The default is to delete missing data }
  \item{skew}{ Should the skew be calculated? }
  \item{ranges}{ Should the range be calculated. }
}
\details{In basic data analysis it is vital to get basic descriptive statistics. 
Procedures such as \code{\link{summary}} and hmisc::describe do so.  The describe function in the \code{\link{psych}} package is meant to produce the most frequently requested stats in psychometric and psychology studies, and to produce them in an easy to read data.frame. The results from describe can be used in graphics functions (e.g., \code{\link{error.crosses}}).

The range statistics (min, max, range) are most useful for data checking to detect coding errors, and should
be found in early analyses of the data.  

The item \code{\link{skew}} is useful to know but will lead to somewhat slower processing times.

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
   standard error\cr
   
 
}
\author{ 
 
 \url{http://personality-project.org/revelle.html} \cr

Maintainer: William Revelle \email{revelle@northwestern.edu} \cr
}
\seealso{ \code{\link{describe.by}}, \code{\link{pairs.panels}}, \code{\link{read.clipboard}}, \code{\link{error.crosses}} }
\examples{

   
#describe(attitude)   #attitude is taken from R data sets

#           var  n  mean    sd median   mad min max range   se
#rating       1 30 64.63 12.17   65.5 10.38  40  85    45 2.22
#complaints   2 30 66.60 13.31   65.0 14.83  37  90    53 2.43
#privileges   3 30 53.13 12.24   51.5 10.38  30  83    53 2.23
#learning     4 30 56.37 11.74   56.5 14.83  34  75    41 2.14
#raises       5 30 64.63 10.40   63.5 11.12  43  88    45 1.90
#critical     6 30 74.77  9.89   77.5  7.41  49  92    43 1.81
#advance      7 30 42.93 10.29   41.0  8.90  25  72    47 1.88

#describe(attitude,skew=TRUE,range=TRUE)

#           var  n  mean    sd median   mad min max range    skew   se
#rating       1 30 64.63 12.17   65.5 10.38  40  85    45 -645.56 2.22
#complaints   2 30 66.60 13.31   65.0 14.83  37  90    53 -508.49 2.43
#privileges   3 30 53.13 12.24   51.5 10.38  30  83    53  694.44 2.23
#learning     4 30 56.37 11.74   56.5 14.83  34  75    41  -87.36 2.14
#raises       5 30 64.63 10.40   63.5 11.12  43  88    45  222.03 1.90
#critical     6 30 74.77  9.89   77.5  7.41  49  92    43 -838.77 1.81
#advance      7 30 42.93 10.29   41.0  8.90  25  72    47  926.20 1.88


## The function is currently defined as
function (x, digits = 2,na.rm=TRUE,skew=FALSE,range=TRUE)
#basic stats after dropping non-numeric data
#somewhat faster if we don't do skews
{    #first, define a local function
    valid <- function(x) {  sum(!is.na(x))}
   		 
    if (is.vector(x) )          #do it for vectors or 
    	{
    	    stats = matrix(rep(NA,6),ncol=6)    #create a temporary array
    
			stats[1, 1] <-  mean(x, na.rm=na.rm )
			stats[1, 2] <-  median(x,na.rm=na.rm  )
			stats[1, 3] <-  min(x, na.rm=na.rm )
			stats[1, 4] <-  max(x, na.rm=na.rm )
			stats[1, 5] <-  skew(x,na.rm=na.rm  )
			stats[1, 6] <-  valid(x )
        
    	}
    	
    	
    else  {len = dim(x)[2]     #do it for matrices or data.frames
    
    stats = matrix(rep(NA,len*6),ncol=6)    #create a temporary array
    for (i in 1:len) {
    	if (is.numeric(x[,i])) {   #just do this for numeric data
			stats[i, 1] <-  mean(x[,i], na.rm=na.rm )
		if (range) {
			stats[i, 2] <-  median(x[,i],na.rm=na.rm  )	
				stats[i, 3] <-  min(x[,i], na.rm=na.rm )
				stats[i, 4] <-  max(x[,i], na.rm=na.rm )}
		    if (skew) {stats[i, 5] <-  skew(x[,i],na.rm=na.rm  )}
			stats[i, 6] <-  valid(x[,i] )
        		}
    	}
    	}
    if (range)
    	{if(skew){temp <-  data.frame(var = seq(1:len),n = stats[,6],mean=stats[,1],
    	sd = sd(x,na.rm=TRUE), median = stats[, 
        2],min= stats[,3],max=stats[,4], range=stats[,4]-stats[,3],skew = stats[, 5])}
         
      	 else {temp <-  data.frame(var = seq(1:len),n = stats[,6],mean=stats[,1],
      	 sd = sd(x,na.rm=TRUE), median = stats[, 
        2],min= stats[,3],max=stats[,4], range=stats[,4]-stats[,3])}}
        
        else {if(skew){temp <-  data.frame(var = seq(1:len),n = stats[,6],mean=stats[,1], 
        sd = sd(x,na.rm=TRUE),skew = stats[, 5])}
       else {temp <-  data.frame(var = seq(1:len),n = stats[,6],mean=stats[,1],
       sd = sd(x,na.rm=TRUE))}}
                
    answer <-  round(data.frame(temp, se = temp$sd/sqrt(temp$n)),  digits)
     return(answer)
  }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
