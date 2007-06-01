\name{principal}
\alias{principal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Principal components analysis }
\description{Does an eigen value decomposition and returns eigen values, loadings, and degree of fit for a specified number of components.  Basically just is doing a principal components for n principal components.  Can show the residual correlations as well. The quality of reduction in the squared correlations is reported by comparing residual correlations to original correlations. Unlike princomp, this returns a subset of just the best nfactors. The eigen vectors are rescaled by the sqrt of the eigen values to produce the component loadings more typical in factor analysis.
}
\usage{
principal(r, nfactors = 0, residuals = FALSE,rotate=FALSE, digits=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{a correlation matrix }
  \item{nfactors}{Number of components to extract }
  \item{residuals}{ FALSE, do not show residuals, TRUE, report residuals }
  \item{rotate} {if TRUE, do varimax rotation, else report unrotated results}
  \item{digits }{digits =2    Accuracy of answers as well as display}
}
\details{Useful for those cases where the correlation matrix is improper (perhaps because of SAPA techniques).
}
\value{
  
  \item{values }{Eigen Values of all components -- useful for a scree plot}
  \item{loadings }{A standard loading matrix}
  \item{fit }{Fit of the model to the correlation matrix }
  \item{residual }{Residual matrix -- if requested}
}

\author{ William Revelle}

\seealso{\code{\link{VSS}},\code{\link{factor2cluster}},\code{\link{factor.pa}}, \code{\link{factor.congruence}}}
\examples{
#Four principal components of the Harmon 24 variable problem
#compare to a four factor principal axes solution using factor.congruence
pc <- principal(Harman74.cor$cov,4,rotate=TRUE)
pa <- factor.pa(Harman74.cor$cov,4,rotate=TRUE)
round(factor.congruence(pc,pa),2)


}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
