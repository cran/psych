\name{factor.fit}
\alias{factor.fit}
\title{   How well does the factor model fit a correlation matrix. Part of the VSS package }
\description{The basic factor or principal components model is that a correlation or covariance matrix may be reproduced by the product of a factor loading matrix times its transpose: F'F or P'P.  One simple index of fit is the 1 - sum squared residuals/sum squared original correlations. This fit index is used by \code{\link{VSS}}, \code{\link{ICLUST}}, etc. 
}
\usage{
factor.fit(r, f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{a correlation matrix }
  \item{f}{A factor matrix of loadings.}
}
\details{There are probably as many fit indices as there are psychometricians.  This fit is a plausible estimate of the amount of reduction in a correlation matrix given a factor model.  Note that it is sensitive to the size of the original correlations.  That is, if the residuals are small but the original correlations are small, that is a bad fit.  
}
\value{fit
}

\author{ William Revelle}


\seealso{ \code{\link{VSS}}, \code{\link{ICLUST}} }
\examples{
\dontrun{
#compare the fit of 4 to 3 factors for the Harman 24 variables
fa4 <- factanal(x,4,covmat=Harman74.cor$cov)
round(factor.fit(Harman74.cor$cov,fa4$loading),2)
#[1] 0.9
fa3 <- factanal(x,3,covmat=Harman74.cor$cov)
round(factor.fit(Harman74.cor$cov,fa3$loading),2)
#[1] 0.88

}

}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
