\name{p.rep}
\alias{p.rep}
\alias{p.rep.f}
\alias{p.rep.t}
\alias{p.rep.r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the probability of replication for an F, t, or r and estimate effect size }
\description{The probability of replication of an experimental or correlational finding as discussed by Peter Killeen (2005) is the probability of finding an effect in the same direction upon an exact replication.  For articles submitted to Psychological Science, p.rep needs to be reported. 

F, t, p and r are all estimates of the size of an effect.  But F, t, and p also are also a function of the sample size.  Effect size, d prime, may be expressed as differences between means compared to within cell standard deviations, or as a correlation coefficient.  These functions convert p, F, and t to d prime and the r equivalent.
}
\usage{
p.rep(p = 0.05, n=NULL,twotailed = FALSE)
p.rep.f(F,df2,twotailed=FALSE) 
p.rep.r(r,n,twotailed=TRUE) 
p.rep.t(t,df,df2=NULL,twotailed=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{conventional probability of statistic (e.g., of F, t, or r)}
  \item{F}{The F statistic}
  \item{df}{Degrees of freedom of the t-test, or of the first group if unequal sizes}
  \item{df2}{Degrees of freedom of the denominator of F or the second group in an unequal sizes t test}
  \item{r}{Correlation coefficient}
  \item{n}{Total sample size if using r }
  \item{t}{t-statistic if doing a t-test or testing significance of a regression slope}
  \item{twotailed}{Should a one or two tailed test be used? }
}
\details{The conventional Null Hypothesis Significance Test (NHST) is the likelihood of observing the data given the null hypothesis of no effect.  But this tells us nothing about the probability of the null hypothesis.  Peter Killeen (2005) introduced the probability of replication as a more useful measure.  The probability of replication is the probability that an exact replication study will find a result in the \emph{same direction} as the original result.

p.rep is based upon a 1 tailed probability value of the observed statistic.  

Other frequently called for statistics are estimates of the effect size, expressed either as Cohen's d, Hedges g, or the equivalent value of the correlation, r. 

For p.rep.t, if the cell sizes are unequal, the effect size estimates are adjusted by the ratio of the mean cell size to the harmonic mean cell size (see Rownow et al., 2000).   
}
\value{
  \item{p.rep }{Probability of replication}
  \item{dprime }{Effect size (Cohen`s d) if more than just p is specified}
  \item{prob}{Probability of F, t, or r.  Note that this can be either the one-tailed or two tailed probability value.}
  \item{r.equivalent}{For t-tests, the r equivalent to the t (see Rosenthal and Rubin(2003), Rosnow, Rosenthal, and Rubin, 2000))}.
}
\references{
Cummings, Geoff (2005) Understanding the average probability of replication: comment on Killeen 2005). Psychological Science, 16, 12, 1002-1004). \cr

Killeen, Peter H. (2005) An alternative to Null-Hypothesis Significance Tests.  Psychological Science, 16, 345-353 \cr

Rosenthal, R. and Rubin, Donald B.(2003), r-sub(equivalent): A Simple Effect Size Indicator.  Psychological Methods, 8, 492-496.


Rosnow, Ralph L., Rosenthal, Robert and Rubin, Donald B. (2000) Contrasts and correlations in effect-size estimation, Psychological Science, 11.  446-453.

}


\note{ The p.rep value is the one tailed probability value of obtaining a result in the same direction.
}

\examples{

p.rep(.05)  #probability of replicating a result if the original study had a  p = .05
p.rep.f(9.0,98)  #probability of replicating  a result with F = 9.0 with 98 df
p.rep.r(.4,50)    #probability of replicating a result if r =.4 with n = 50
p.rep.t(3,98)   #probability of replicating a result if t = 3 with df =98
p.rep.t(2.14,84,14) #effect of equal sample sizes (see Rosnow et al.)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{univar}

