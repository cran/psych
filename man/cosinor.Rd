\name{cosinor}
\alias{cosinor}
\alias{circadian.mean}
\alias{circular.mean}
\alias{circadian.cor}
\alias{circular.cor} 
\alias{circadian.linear.cor}
\title{Functions for analysis of circadian or diurnal data}
\description{Circadian data are periodic with a phase of 24 hours. These functions  find the best fitting phase angle (cosinor), the circular mean,  circular correlation with circadian data, and the linear by circular correlation} 

\usage{
cosinor(angle,x=NULL,code=NULL,period=24,plot=FALSE,opti=FALSE)
circadian.mean(angle, hours=TRUE)
circular.mean(angle) #angles in radians
circadian.cor(angle,hours=TRUE)  #angles in radians
circular.cor(angle) #angles in radians
circadian.linear.cor(angle,x,hours=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{angle}{A data frame or matrix of observed values with the time of day as the first value (unless specified in code) angle can be specified either as hours or as radians)}
  \item{code}{A subject identification variable}
  \item{period}{Although time of day is assumed to have a 24 hour rhythm, other rhythms may be fit. }
 \item{plot}{if TRUE, then plot the first variable (angle)} 
 \item{opti}{opti=TRUE: iterative optimization (slow) or opti=FALSE: linear fitting (fast)}
 \item{hours}{If TRUE, measures are in 24 hours to the day, otherwise, radians}
 \item{x}{A set of external variables to correlate with the phase angles}
}
\details{
When data represent angles (such as the hours of peak alertness or peak tension during the day), we need to apply circular statistics rather than the more normal linear statistics (see Jammalamadaka (2006) for a very clear set of examples of circular statistics). The generalization of the mean to circular data is to convert each angle into a vector, average the x and y coordinates, and convert the result back to an angle. The generalization of Pearson correlation to circular statistics is straight forward and is implemented in cor.circular in the circular package and in \code{\link{circadian.cor}} here.  Just as the Pearson r is a ratio of covariance to the square root of the product of two variances, so is the circular correlation.  The circular covariance of two circular vectors is defined as the average product of the sines of the deviations from the circular mean.  The variance is thus the average squared sine of the angular deviations from the circular mean.  Circular statistics are used for data that vary over a period (e.g., one day) or over directions (e.g., wind direction or bird flight).  Jammalamadaka and Lund (2006)  give a very good example of the use of circular statistics in calculating wind speed and direction.  

The code from CircStats and circular was adapted to allow for analysis of data from various studies of mood over the day.

\code{\link{circular.mean}} and \code{\link{circular.cor}} are just \code{\link{circadian.mean}} and \code{\link{circadian.cor}}but with input given in radians rather than hours.

The cosinor function will either iteratively fit cosines of the angle to the observed data (opti=TRUE) or use the circular by linear regression to estimate the best fitting phase angle.  If cos.t <- cos(time) and sin.t = sin(time) (expressed in hours), then beta.c and beta.s may be found by regression and the phase is \eqn{sign(beta.c) * acos(beta.c/\sqrt(beta.c^2 + beta.s^2)) * 12/pi}

Simulations (see examples) suggest that with incomplete times, perhaps the optimization procedure yields slightly better fits with the correct phase than does the linear model, but the differences are very small. In the presence of noisey data, these advantages seem to reverse.  The recommendation thus seems to be to use the linear model approach (the default).
}
\value{
  \item{phase }{The phase angle that best fits the data}
  \item{fit }{Value of the correlation of the fit}
  \item{mean.angle}{A vector of mean angles}
  \item{R}{A matrix of circular correlations or linear by circular correlations}

}
\references{ See circular statistics 
Jammalamadaka, Sreenivasa and Lund, Ulric (2006),The effect of wind direction on ozone levels: a case study, Environmental and Ecological Statistics, 13, 287-298.	
}
\author{William Revelle }

\seealso{See the circular and CircStats packages. }

\examples{
time <- seq(1:24)
pure <- matrix(time,24,18)
pure <- cos((pure + col(pure))*pi/12)
matplot(pure,type="l",main="Pure circadian arousal rhythms",xlab="time of day",ylab="Arousal") 
p <- cosinor(time,pure)
#set.seed(42)
noisy <- pure + rnorm(24*18)
n <- cosinor(time,noisy) 
#small.pure <- pure[c(6:18),]
small.pure <- pure[c(8,11,14,17,20,23),]
#small.noisy <- noisy[c(6:18),]
small.noisy <- noisy[c(8,11,14,17,20,23),]
matplot(small.noisy,type="l",main="Noisy circadian arousal rhythms",
         xlab="time of day",ylab="Arousal") 
#sp <- cosinor(time[c(6:18)],small.pure)  #linear fit
sp <- cosinor(time[c(8,11,14,17,20,23)],small.pure)
spo <- cosinor(time[c(8,11,14,17,20,23)],small.pure,opti=TRUE) #iterative fit
sn <- cosinor(time[c(8,11,14,17,20,23)],small.noisy) #linear
sno <- cosinor(time[c(8,11,14,17,20,23)],small.noisy,opti=TRUE) #iterative
sum.df <- data.frame(pure=p,noisy = n, small=sp,small.noise = sn, 
        small.opt=spo,small.noise.opt=sno)
round(sum.df,2)
round(circadian.cor(sum.df[,c(1,3,5,7,9,11)]),2)  #compare alternatives 
round(cor(sum.df[,c(2,4,6,8,10,12)]),2) 
}

\keyword{ multivariate }
