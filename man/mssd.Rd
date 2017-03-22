\name{mssd}
\alias{mssd}
\alias{rmssd}
\alias{autoR}
\title{Find von Neuman's Mean Square of Successive Differences}

\description{
Von Neuman et al. (1941) discussed the Mean Square of Successive Differences as a measure of variability that takes into account gradual shifts in mean. This is appropriate when studying errors in ballistics or variability and stability in mood when studying affect. For random data, this will be twice the variance, but for data with a sequential order and a positive autocorrelation, this will be much smaller. Since the mssd is just twice the variance - the autocorrelation, it is thus possible to also find the autocorrelation for a particular lag. 
}
\usage{
mssd(x,group=NULL, lag = 1,na.rm=TRUE)
rmssd(x,group=NULL, lag=1, na.rm=TRUE)
autoR(x,group=NULL,lag=1,na.rm=TRUE)
}

\arguments{
  \item{x}{a vector, data.frame or matrix}
  \item{lag}{the lag to use when finding \code{\link{diff}}     }
  \item{group}{A column of the x data.frame to be used for grouping}
  \item{na.rm}{Should missing data be removed?}
}
\details{When examining multiple measures within subjects, it is sometimes useful to consider the variability of trial by trial observations in addition to the over all between trial variation.  The Mean Square of Successive Differences (mssd) and root mean square of successive differences (rmssd) find the variance or standard deviation of the trial to trial differences. 

\eqn{\sigma^2 = \Sigma(x_i - x_{i+1})^2 /(n-1) } 

In the case of multiple subjects  (groups) with multiple observations per subject (group), specify the grouping variable will produce output for each group.   


Similar functions are available in the matrixStats package. This is just the variance and standard deviation applied to the result from the \code{\link{diff}} function.

Perhaps useful when studying mood, the \code{\link{autoR}} function finds the autocorrelation for each item for the specified lag.  It also returns the rmssd (root means square successive difference).  


}
\value{
The average squared successive difference (mssd) and the square root of the average squared successive difference (rmssd).  Note that this is not the same as the standard deviation of the lagged differences.
}
\references{

Jahng, Seungmin and Wood, Phillip K and Trull, Timothy J. Analysis of affective instability in ecological momentary assessment: Indices using successive difference and group comparison via multilevel modeling. Psychological methods (2008) 13, 354-375. 

Von Neumann, J., Kent, R., Bellinson, H., and Hart, B. (1941). The mean square successive difference. The Annals of Mathematical Statistics, pages 153-162.


}
\author{William Revelle}


\seealso{See Also \code{\link{rmssd}} for the standard deviation or  \code{\link{describe}} for more conventional statistics.  \code{\link{describeBy}} and \code{\link{statsBy}} give group level statistics.  See also \code{link{mlr}}, \code{link{mlreliability}}, \code{link{mlPlot}} for other ways of examining within person variability over time. 
}
\examples{
t <- seq(-pi, pi, .1)
trial <- 1: length(t)
gr <- trial \%\% 8 
c <- cos(t)
ts <- sample(t,length(t))
cs <- cos(ts)
x.df <- data.frame(trial,gr,t,c,ts,cs)
rmssd(x.df)
rmssd(x.df,gr)
autoR(x.df,gr)
describe(x.df)
#pairs.panels(x.df)
#mlPlot(x.df,grp="gr",Time="t",items=c(4:6))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}
