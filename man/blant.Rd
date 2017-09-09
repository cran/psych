\name{blant}
\alias{blant}
\docType{data}
\title{A 29 x 29 matrix that produces weird factor analytic results}
\description{Normally, min.res factor analysis and maximum likelihood produce very similar results.  This data set (from Alexandra Blant) does not.  Warnings are given for the min.res solution, the pa solution, but not the old.min nor the mle solution.  Included as a test case for the factor analysis function.
}
\usage{data("blant")}
\format{
  The format is:
 num [1:29, 1:29] 1 0.77 0.813 0.68 0.717 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:29] "V1" "V2" "V3" "V4" ...
}
\details{
This data matrix was sent by Alexandra Blant as an example of a problem with the minres solution in the \code{\link{fa}} function.  The default solution, using fm="minres" issues a warning that the solution has improper factor score weights.  This is not the case for the fm="old.min" and fm="mle" options, but is for fm="pa", fm="ols".

The residuals are indeed smaller for fm="minres" than for fm="old.min" or fm="mle".

"old.min" attempts to find the minimum residual but uses the gradient for mle. This was the approach until version 1.7.5 but was changed (see the help page for fa) following extensive communication with Hao Wu.

The problem with this matrix is probably that it is almost singular, with some smcs approaching 1 and the smallest three eigenvalues of .006, .004 and .001.

This problem matrix was provided by Alexandra Blant.  
}
\source{Alexandra Blant, personal communication}


\examples{
data(blant)
#compare
f5 <- fa(blant,5,rotate="none")  #the default minres 
f5.old <-fa(blant,5, fm="old.min",rotate="none") #old version of minres
f5.mle <- fa(blant,5,fm="mle",rotate= "none") #maximum likelihood
#compare solutions
factor.congruence(list(f5,f5.old,f5.mle))
#compare sums of squared residuals
sum(residuals(f5,diag=FALSE)^2,na.rm=TRUE) # 1.355489
sum(residuals(f5.old,diag=FALSE)^2,na.rm=TRUE) # 1.539757
sum(residuals(f5.mle,diag=FALSE)^2,na.rm=TRUE) # 2.402092

#but, when we divide the squared residuals by the original (squared) correlations, we find 
#a different ordering of fit
f5$fit     # 0.9748177
f5.old$fit  # 0.9752774
f5.mle$fit   # 0.9603324

}
\keyword{datasets}
