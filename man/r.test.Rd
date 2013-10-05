\name{r.test}
\alias{r.test}
\title{Tests of significance for correlations}
\description{Tests the significance of a single correlation, the difference between two independent correlations, the difference between two dependent correlations sharing one variable (Williams's Test), or the difference between two dependent correlations with different variables (Steiger Tests).
}
\usage{
r.test(n, r12, r34 = NULL, r23 = NULL, r13 = NULL, r14 = NULL, r24 = NULL, 
          n2 = NULL,pooled=TRUE, twotailed = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size of first group }
  \item{r12}{Correlation to be tested}
  \item{r34}{Test if this correlation is different from r12, if r23 is specified, but r13 is not, then r34 becomes r13  }
  \item{r23}{if ra = r(12) and rb = r(13)  then test for differences of dependent correlations given r23}
  \item{r13}{implies  ra =r(12) and rb =r(34)  test for difference of dependent correlations }
  \item{r14}{implies   ra =r(12) and rb =r(34) }
  \item{r24}{ ra =r(12) and rb =r(34)}
  \item{n2}{n2 is specified in the case of two independent correlations. n2 defaults to n if if not specified }
  \item{pooled}{use pooled estimates of correlations}
  \item{twotailed}{ should a twotailed or one tailed test be used }
}
\details{Depending upon the input, one of four different tests of correlations is done.
1) For a sample size n, find the t value for a single correlation. 

2) For sample sizes of n and n2 (n2 = n if not specified) find the z of the difference between the z transformed correlations divided by the standard error of the difference of two z scores. 

3) For sample size n, and correlations r12, r13 and r23 test for the difference of two dependent correlations (r12 vs r13).  

4) For sample size n, test for the difference between two dependent correlations involving different variables.

For clarity, correlations may be specified by value.  If specified by location and if doing the test of dependent correlations, if three correlations are specified, they are assumed to be in the order r12, r13, r23. Consider the example  the example from Steiger: 
where Masculinity at time 1 (M1) correlates with Verbal Ability .5 (r12), femininity at time 1 (F1) correlates with Verbal ability  r13 =.4, and M1 correlates with F1 (r23= .1).  Then, given the correlations: r12 = .4, r13 = .5, and r23 = .1, t = -.89 for n =103, i.e.,
r.test(n=103, r12=.4, r13=.5,r23=.1)  



}
\value{
  
  \item{test}{Label of test done}
  \item{z}{z value for tests 2 or 4}
  \item{t}{t value for tests 1 and 3}
  \item{p}{probability value of z or t}

}
\references{Olkin, I. and Finn, J. D. (1995). Correlations redux. Psychological Bulletin, 118(1):155-164. 

Steiger, J.H. (1980), Tests for comparing elements of a correlation matrix, Psychological Bulletin, 87, 245-251. 

Williams, E.J.  (1959) Regression analysis. Wiley, New York, 1959.


}
\author{William Revelle }
\note{Steiger specifically rejects using the Hotelling T test to test the difference between correlated correlations. Instead, he recommends Williams' test. (See also Dunn and Clark, 1971). These tests follow Steiger's advice.

}
\seealso{ See also  \code{\link{corr.test}} which tests all the elements of a  correlation matrix, and \code{\link{cortest.mat}} to compare two matrices of correlations.  r.test  extends the tests in \code{\link{paired.r}},\code{\link{r.con}}}
\examples{

n <- 30 
r <- seq(0,.9,.1) 
rc <- matrix(r.con(r,n),ncol=2) 
test <- r.test(n,r)
r.rc <- data.frame(r=r,z=fisherz(r),lower=rc[,1],upper=rc[,2],t=test$t,p=test$p) 
round(r.rc,2) 

r.test(50,r)
r.test(30,.4,.6)       #test the difference between two independent correlations
r.test(103,.4,.5,.1)   #Steiger case A of dependent correlations 
r.test(n=103, r12=.4, r13=.5,r23=.1)  
#for complicated tests, it is probably better to specify correlations by name
r.test(n=103,r12=.5,r34=.6,r13=.7,r23=.5,r14=.5,r24=.8)   #steiger Case B 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
