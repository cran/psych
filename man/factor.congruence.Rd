\name{factor.congruence}
\alias{factor.congruence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coefficient of factor congruence }
\description{Given two sets of factor loadings, report their degree of congruence. 
}
\usage{
factor.congruence(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix of factor loadings }
  \item{y}{ A second matrix of factor loadings }

}
\details{Find the coefficient of factor congruence between two sets of factor loadings. 

It is an interesting exercise to compare factor congruences with the correlations of factor loadings.  Factor congruences are based upon the raw cross products, while correlations are based upon centered cross products.  

Input may either be matrices or factor analysis output (which includes a loadings object), or a mixture of the two.

Factor congruences are the cosines of pairs vectors defined by the loadings matrix and based at the origin.  The correlations of factor loadings are cinoses of the vectors based at the mean loading for each factor.  

}
\value{A matrix of factor congruences.
}

\references{Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates.
\cr
Revelle, W. (In preparation) An Introduction to Psychometric Theory with applications in R (\url{http://personality-project.org/r/book/})
}
\author{
 \email{revelle@northwestern.edu} \cr
 \url{http://personality-project.org/revelle.html}}

\seealso{ \code{\link{principal}}, \code{\link{factor.pa}}}
 
\examples{
#fa <- factanal(x,4,covmat=Harman74.cor$cov)
#pc <- principal(Harman74.cor$cov,4)
#pcv <- varimax(pc$loading)
#round(factor.congruence(fa,pcv),2)
#
#round(factor.congruence(pcv,fa),2)
#    Factor1 Factor2 Factor3 Factor4
#PC1    1.00    0.60    0.45    0.55
#PC2    0.44    0.49    1.00    0.56
#PC3    0.54    0.99    0.44    0.55
#PC4    0.47    0.52    0.48    0.99


#compare with 
#round(cor(fa$loading,pcv$loading),2)


}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
