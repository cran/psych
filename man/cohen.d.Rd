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
\alias{m2d}
\alias{d2OVL}
\alias{d2OVL2}
\alias{d2CL}
\alias{d2U3}
\alias{cd.validity}


\title{Find Cohen d and confidence intervals} 

\description{
Given a data.frame or matrix, find the standardized mean difference (Cohen's d) and confidence intervals for each variable depending upon a grouping variable.  Convert the d statistic to the r equivalent, report the student's t statistic and associated p values, and return statistics for both values of the grouping variable.  The Mahalanobis distance between the centroids of the two groups in the space defined by all the variables ia also found.  Confidence intervals for Cohen d for one group (difference from 0) may also be found. Several measures of the distributional overlap (e.g. OVL, OVL2, etc.) are available.
}
\usage{
cohen.d(x, group,alpha=.05,std=TRUE,sort=NULL,dictionary=NULL,MD=TRUE,data=NULL)
d.robust(x,group,trim=.2)
cohen.d.ci(d,n=NULL,n2=NULL,n1=NULL,alpha=.05)
d.ci(d,n=NULL,n2=NULL,n1=NULL,alpha=.05)
cohen.d.by(x,group,group2,alpha=.05,MD=TRUE)
d2r(d)
r2d(rho)
d2t(d,n=NULL,n2=NULL,n1=NULL)
t2d(t,n=NULL,n2=NULL,n1=NULL)
m2t(m1,m2,s1,s2,n1=NULL,n2=NULL,n=NULL,pooled=TRUE)  #returns d invisibily
m2d(m1,m2,s1,s2,n1=NULL,n2=NULL,n=NULL,pooled=TRUE)
d2OVL(d)  #Percent overlap for 1 distribtion
d2OVL2(d)  #Percent overlap joint distribution
d2CL(d)   #Common language effect size
d2U3(d)   #Proportion in higher group exceedding median of lower group
cd.validity(d, keys, abs=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data frame or matrix  (can be specified in formula mode) }
  \item{group}{Some dichotomous grouping variable (may be specified using formula input (see example))}
  \item{group2}{Apply cohen.d for each of the subgroups defined by group2 (may be specified by formula as well)}
  \item{data}{If using formula mode and specifying a particular variable (see example)}
 \item{d}{An effect size}
 \item{keys}{A list of scoring keys (similar to scoreItems)}
 \item{abs}{When finding average cd validities, should we take absolute values (TRUE)}
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
 \item{MD}{Find Mahalanobis distance in cohen.d.}
}
\details{
There are many ways of reporting how two groups differ.  Cohen's d statistic is just the differences of means expressed in terms of the pooled within group standard deviation.  This is insensitive to sample size.  r is the a universal measure of effect size that is a simple function of d, but is bounded -1 to 1.  The t statistic is merely d * sqrt(n)/2 and thus reflects sample size.   

\deqn{d  = \frac{M2 - M1}{Sp}}{(M2- M1)/Sp}  
where Sp is the pooled standard deviation.  

\deqn{Sp  = \sqrt{\frac{(n1-1)*s1^2 + (n2-1)* s2^2}{N}  }  }{\sqrt{((n1-1)*s1^2 + (n2-1)* s2^2)/{N}  }  }

Cohens d uses N as the divisor for the pooled sums of squares.  Hedges g uses N-2.  

Confidence intervals for Cohen's d may be found by converting the d to a t, finding the confidence intervals for t, and then converting those back to ds.  This take advantage of the uniroot function and the non-centrality parameter of the t distribution.


The results of \code{\link{cohen.d}} may be displayed using the \code{\link{error.dots}} function.  This will include the labels provided in the dictionary.  

In the case of finding the confidence interval (using \code{\link{cohen.d.ci}} for a comparison against 0 (the one sample case), specify n1.  This will yield a d = t/sqrt(n1)  whereas in the case of the difference between two samples, d = 2*t/sqrt(n) (for equal sample sizes n = n1+ n2) or d = t/sqrt(1/n1 + 1/n2)  for the case of unequal sample sizes.

Since we find d and then convert this to t, using d2t,  the question is how to pool the variances. Until 7/14/21 I was using the total n to estimate the t and thus the p values.  In response to a query (see news), I switched to using the actual sample size ns (n1 and n2) and then finding t based upon the hedges g value.  This produces t values as reported by t.test with the var.equal = TRUE option.

It is probably useful to comment that the various confidence intervals reported are based upon normal theory and should be interpreted cautiously.  


\code{\link{cohen.d.by}} will find Cohen's d for groups for each subset of the data defined by group2.  The summary of the output produces a simplified listing of the d values for each variable for each group.  May be called directly from cohen.d by using formula input and specifying two grouping variables. 

\code{\link{d.robust}} follows Algina et al. 2005) to find trimmed means (trim =.2) and Winsorize variances (trim =.2).  Supposedly, this provides a more robust estimate of effect sizes.

\code{\link{m2t}} reports Student's t.test for two groups given their means, standard deviations, and sample size.  This is convenient when checking statistics where those estimates are provided, but the raw data are not available.  By default, it gives the pooled estimate of variance, but if pooled is FALSE, it applies Welch's correction.

The Mahalanobis Distance combines the individual ds and weight them by their unique contribution:  \eqn{D = \sqrt{d' R^{-1}d}}{D = \sqrt{d' R^{-1}d}}.
 By default, \code{\link{cohen.d}} will find the Mahalanobis distance between the two groups (if there is more than one DV.)  This requires finding the correlation of all of the DVs and can fail if that matrix is not invertible because some pairs do not exist.  Thus, setting MD=FALSE will prevent the Mahalanobis calculation.
 
 Marco del Giudice (2019) has a very helpful paper discussing how to interpret d and Md in terms of various overlap coefficients. These may be found by the use of the \code{\link{d2OVL}} (percent overlap for 1 distribution), \code{\link{d2OVL2}} percent overlap of joint distributions, \code{\link{d2CL}} (the common language effect size), and \code{\link{d2U3}} (proportion in higher group exceeding median of the lower group).
 
 \deqn{OVL  = 2\phi(-d/2)}{OVL = 2 \phi(-d/2)}  is the proportion of overlap (and gets smaller the larger the d).
 where  Phi  is the cumulative density function of the normal distribution.
 \deqn{OVL_2  = \frac{OVL}{2-OVL}}{OVL_2 = OVL/(2-OVL)}

The proportion of individuals in one group above the median of the other group is U3
\deqn{U_3 = \phi_d}{U_3 = \phi(d)}.

The Common Language Effect size 
 \deqn{CL  = \phi (d *\sqrt{2})}{CL = \phi(d * \sqrt(2) )}
 These last two get larger with (abs (d)).  
 For graphic displays of Cohen's d and Mahalanobis D, see the \code{\link{scatterHist}} examples, or the example from the psychTools::GERAS data set. 

}
\value{
\item{d}{Cohen's d statistic, including the upper and lower confidence levels}
\item{hedges.g}{Hedge's g statistic, including the upper and lower confidence levels}
\item{M.dist}{Mahalanobis distance between the two groups}
\item{t}{Student's t statistic}
\item{r}{The point biserial r equivalent of d}
\item{n}{sample size used for each analysis}
\item{p}{The probability of abs(t)>0}
\item{descriptive}{The descriptive statistics for each group}
\item{OVL}{etc. some of the measures of overlap discussed by DelGiudice, 2009}

}
\note{Cohen and Hedges differ in they way they calculate the pooled within group standard deviation. I find the treatment by McGrath and Meyer to be most helpful in understanding the differences.
}
\references{Cohen, Jackob (1988) Statistical Power Analysis for the Behavioral Sciences. 2nd Edition, Lawrence Erlbaum Associates. 

Algina, James and Keselman, H. J. and Penfield, Randall D. (2005) An Alternative to Cohen's Standardized Mean Difference Effect Size: A Robust Parameter and Confidence Interval in the Two Independent Groups Case. Psychological Methods.  10, 317-328.

Goulet-Pelletier, Jean-Christophe and Cousineau, Denis.(2018) A review of effect sizes and their confidence intervals, Part I: The Cohen's d family. The Quantitative Methods for Psychology, 14, 242-265.

Marco Del Giudice (2019) Measuring Sex Differences and Similarities, (in VanderLaan and Wong (ed. ) Gender and sexuality development: Contemporary theory and research.)

McGrath, Robert E and Meyer, Gregory J. (2006) When effect sizes disagree: the case of r and d. Psychological methods, 11, 4, 386-401.



}

\author{
William Revelle
}



\seealso{
\code{\link{describeBy}}, \code{\link{describe}} \code{\link{error.dots}} to display the results.  \code{\link{scatterHist}} to show d and MD for pairs of variables. (See in particular the use of \code{\link{scatterHist}} on psychTools::GERAS daa set.)
}
\examples{
cohen.d(sat.act,"gender")
#robust version
round(d.robust(sat.act,"gender")$robust.d,2)

#formula input is nicer
cohen.d(sat.act ~ gender) #formula input version

#report cohen.d by another group
cd <- cohen.d.by(sat.act,"gender","education")
cohen.d(SATV + SATQ ~ gender, data=sat.act) #just choose two variables
summary(cd)  #summarize the output

#formula version combines these functions
cd <- cohen.d(sat.act ~ gender + education)  #find d by gender for each level of education
summary(cd)

#now show several examples of confidence intervals
#one group (d vs 0)
#consider the t from the cushny data set
t2d( -4.0621,n1=10)
d.ci(-1.284549,n1=10)  #the confidence interval of the effect of drug on sleep
#two groups
d.ci(.62,n=64)  #equal group size
d.ci(.62,n1=35,n2=29) #unequal group size
#several examples of d and t from data
m2d(52.58,-70.65,49.9,47.5) #Terman and Miles 1936

#graphically show the various overlap statistics
curve(d2OVL2(x),0,3,xlab="d",ylab="",lty="dashed",
       main="Four representations of effect size (d) ")
curve(d2OVL(x),0,3,xlab="d",add=TRUE,)
curve(d2CL(x),0,3,add=TRUE)
curve(d2U3(x), add=TRUE,lty="dotted")
text(1,.37,"OVL2")
text(2,.37,"OVL")
text(1,.88,"U3")
text(2, .88,"CL")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
