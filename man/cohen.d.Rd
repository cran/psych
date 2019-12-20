\name{cohen.d}
\alias{cohen.d}
\alias{d.robust}
\alias{cohen.d.ci}
\alias{d.ci}
\alias{cohen.d.by}
\alias{d2r}
\alias{r2d}
\alias{d2t}
\alias{t2d}
\alias{m2t}

\title{Find Cohen d and confidence intervals} 

\description{
Given a data.frame or matrix, find the standardized mean difference (Cohen's d) and confidence intervals for each variable depending upon a grouping variable.  Convert the d statistic to the r equivalent, report the student's t statistic and associated p values, and return statistics for both values of the grouping variable.  The Mahalanobis distance between the centroids of the two groups in the space defined by all the variables ia also found.  Confidence intervals for Cohen d for one group (difference from 0) may also be found. 
}
\usage{
cohen.d(x, group,alpha=.05,std=TRUE,sort=NULL,dictionary=NULL)
d.robust(x,group,trim=.2)
cohen.d.ci(d,n=NULL,n2=NULL,n1=NULL,alpha=.05)
d.ci(d,n=NULL,n2=NULL,n1=NULL,alpha=.05)
cohen.d.by(x,group,group2,alpha=.05)
d2r(d)
r2d(rho)
d2t(d,n=NULL,n2=NULL,n1=NULL)
t2d(t,n=NULL,n2=NULL,n1=NULL)
m2t(m1,m2,s1,s2,n1=NULL,n2=NULL,n=NULL,pooled=TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data frame or matrix}
  \item{group}{Some dichotomous grouping variable}
  \item{group2}{Apply cohen.d for each of the subgroups defined by group2}
 \item{d}{An effect size}
 \item{trim}{The amount of trimming used in finding the means and sds in d.robust}
 \item{n}{Total sample size (of groups 1 and 2)}
 \item{n1}{Sample size of group 1 (if only one group)}
 \item{n2}{Sample size of group 2 }
 \item{pooled}{Pool the two variances}
 \item{t}{Student's t statistic}
 \item{alpha}{1-alpha is the width of the confidence interval}
 \item{std}{Find the correlation rather covariance matrix}
 \item{rho}{A correlation to be converted to an effect size}
 \item{m1}{Mean of group 1}
 \item{m2}{Mean of group 2}
 \item{s1}{Standard deviation of group 1}
 \item{s2}{Standard deviation of group 2}
 \item{sort}{Should we sort (and if so, in which direction), the results of cohen.d? 
 Directions are "decreasing" or  "increasing".}
 \item{dictionary}{What are the items being described?}
}
\details{
There are many ways of reporting how two groups differ.  Cohen's d statistic is just the differences of means expressed in terms of the pooled within group standard deviation.  This is insensitive to sample size.  r is the a universal measure of effect size that is a simple function of d, but is bounded -1 to 1.  The t statistic is merely d * sqrt(n)/2 and thus reflects sample size.   

Confidence intervals for Cohen's d may be found by converting the d to a t, finding the confidence intervals for t, and then converting those back to ds.  This take advantage of the uniroot function and the non-centrality parameter of the t distribution.

The results of \code{\link{cohen.d}} may be displayed using the \code{\link{error.dots}} function.  This will include the labels provided in the dictionary.  

In the case of finding the confidence interval for a comparison against 0 (the one sample case), specify n1.  This will yield a d = t/sqrt(n1)  whereas in the case of the differnece between two samples, d = 2*t/sqrt(n) (for equal sample sizes n = n1+ n2) or d = t/sqrt(1/n1 + 1/n2)  for the case of unequal sample sizes.
 
\code{\link{cohen.d.by}} will find Cohen's d for groups for each subset of the data defined by group2.  The summary of the output produces a simplified listing of the d values for each variable for each group.

\code{\link{d.robust}} follows Algina et al. 2005) to find trimmed means (trim =.2) and Winsorize variances (trim =.2).  Supposedly, this provides a more robust estimate of effect sizes.

\code{\link{m2t}} reports Student's t.test for two groups given their means, standard deviations, and sample size.  This is convenient when checking statistics where those estimates are provided, but the raw data are not available.  By default,it gives the pooled estimate of variance, but if pooled is FALSE, it applies Welch's correction.
 
}
\value{
\item{d}{Cohen's d statistic, including the upper and lower confidence levels}
\item{hedges.g}{Hedge's g statistic}
\item{M.dist}{Mahalanobis distance between the two groups}
\item{t}{Student's t statistic}
\item{r}{The point biserial r equivalent of d}
\item{n}{sample size used for each analysis}
\item{p}{The probability of abs(t)>0}
\item{descriptive}{The descriptive statistics for each group}

}
\references{Cohen, Jackob (1988) Statistical Power Analysis for the Behavioral Sciences. 2nd Edition, Lawrence Erlbaum Associates. 

Algina, James and Keselman, H. J. and Penfield, Randall D. (2005) An Alternative to Cohen's Standardized Mean Difference Effect Size: A Robust Parameter and Confidence Interval in the Two Independent Groups Case. Psychological Methods.  10, 317-328.

}
\author{
William Revelle
}



\seealso{
\code{\link{describeBy}}, \code{\link{describe}}
}
\examples{
cohen.d(sat.act,"gender")
cd <- cohen.d.by(sat.act,"gender","education")
summary(cd)  #summarize the output

#now show several examples of confidence intervals
#one group (d vs 0)
#consider the t from the cushny data set
t2d( -4.0621,n1=10)
d.ci(-1.284549,n1=10)  #the confidence interval of the effect of drug on sleep
#two groups
d.ci(.62,n=64)  #equal group size
d.ci(.62,n1=35,n2=29) #unequal group size


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
