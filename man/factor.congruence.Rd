\name{factor.congruence}
\alias{factor.congruence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coefficient of factor congruence }
\description{Given two sets of factor loadings, report their degree of congruence (vector cosine). 
}
\usage{
factor.congruence(x, y,digits=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix of factor loadings }
  \item{y}{ A second matrix of factor loadings }
  \item{digits}{Round off to digits}

}
\details{Find the coefficient of factor congruence between two sets of factor loadings. 

Factor congruences are the cosines of pairs vectors defined by the loadings matrix and based at the origin.  Thus, for loadings that differ only by a scaler (e.g. the size of the eigen value), the factor congruences will be 1.

It is an interesting exercise to compare factor congruences with the correlations of factor loadings.  Factor congruences are based upon the raw cross products, while correlations are based upon centered cross products. That is,  correlations of factor loadings are cosines of the vectors based at the mean loading for each factor.   

Input may either be matrices or factor analysis or principal components analyis output (which includes a loadings object), or a mixture of the two.



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
#pa <- factor.pa(Harman74.cov$cor,4)
#factor.congruence(fa,pc)
#
#    Factor1 Factor2 Factor3 Factor4
#PC1    1.00    0.60    0.45    0.55
#PC2    0.44    0.49    1.00    0.56
#PC3    0.54    0.99    0.44    0.55
#PC4    0.47    0.52    0.48    0.99

# pa <- factor.pa(Harman74.cor$cov,4)
# factor.congruence(fa,pa)
#         PA1  PA3  PA2  PA4
#Factor1 1.00 0.61 0.46 0.55
#Factor2 0.61 1.00 0.50 0.60
#Factor3 0.46 0.50 1.00 0.57
#Factor4 0.56 0.62 0.58 1.00


#compare with 
#round(cor(fa$loading,pc$loading),2)


}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
