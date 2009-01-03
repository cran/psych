\name{cosinor}
\alias{cosinor}
\alias{circadian.mean}
\alias{circadian.cor}
\alias{circadian.linear.cor}
\title{Functions for analysis of circadian or diurnal data}
\description{Circadian data are periodic with a phase of 24 hours. These functions  find the best fitting phase angle (cosinor), the circular mean,  circular correlation with circadian data, and the linear by circular correlation} 

\usage{
cosinor(data, code = NULL, period = 24)
circadian.mean(angle, hours=TRUE)
circadian.cor(angle, hours=TRUE)
circadian.linear.cor(angle,x,hours=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data frame or matrix of observed values with the time of day as the first value (unless specified in code }
  \item{code}{Time of day of measurements}
  \item{period}{Although time of day is assumed to have a 24 hour rhythm, other rhythms may be fit. }
  \item{angle}{A vector, matrix, or data frame of phase angles (either as hours or as radians).}
  \item{hours}{If TRUE, measures are in 24 hours to the day, otherwise, radians}
  \item{x}{A set of external variables to correlate with the phase angles}
}
\details{
When data represent angles (such as the hours of peak alertness or peak tension during the day), we need to apply circular statistics rather than the more normal linear statistics (see Jammalamadaka(2006) for a very clear set of examples of circular statistics). The generalization of the mean to circular data is to convert each angle into a vector, average the x and y coordinates, and convert the result back to an angle. The gneralization of Pearson correlation to circular statistics is straight forward and is implemented in cor.circular in the circular package and in circadian.cor here.  Just as the Pearson r is a ratio of covariance to the square root of the product of two variances, so is the circular correlation.  The circular covariance of two circular vectors is defined as the average product of the sines of the deviations from the circular mean.  The variance is thus the average squared sine of the angular deviations from the circular mean.  Circular statistics are used for data that vary over a period (e.g., one day) or over directions (e.g., wind direction or bird flight).  Jammalamadaka and Lund (2006)  givesa very good example of the use of circular statistics in calculating wind speed and direction.   
The code from CircStats and circular was adapted to allow for an analysis of data from an experimental study of mood over the day.
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{phase }{The phase angle that best fits the data}
  \item{fit }{Value of the correlation of the fit}
  \item{mean.angle}{A vector of mean angles}
  \item{R}{A matrix of circular correlations or linear by circular correlations}

}
\references{ See circular statistics 
Jammalamadaka, Sreenivasa and Lund, Ulric (2006),The effect of wind direction on ozone levels: a case study ,Environmental and Ecological Statistics, 13, 287-298.	
}
\author{William Revelle }

\seealso{See the circular and CircStats packages. }


\keyword{ multivariate }
