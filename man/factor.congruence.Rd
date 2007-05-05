\name{factor.congruence}
\alias{factor.congruence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coefficient of factor congruence }
\description{Given two sets of factor loadings, report their degree of congruence. 
}
\usage{
factor.congruence(x, y, loading=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix of factor loadings }
  \item{y}{ A second matrix of factor loadings }
  \item{loading}{Is the input a loading matrix from a factor analysis or principal components, or is it just a matrix}
}
\details{Find the coefficient of factor congruence between two sets of factor loadings. 

It is an interesting exercise to compare factor congruences with the correlations of factor loadings.  Factor congruences are based upon the raw cross products, while correlations are based upon centered cross products.  

}
\value{A matrix of factor congruences.
}

\references{Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates. Revelle, W. In preparation) An Introduction to Psychometric Theory with applications in R (\url{http://personality-project.org/r/book/})
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



## The function is currently defined as
#function (x,y,loading=TRUE) {
# if (loading) {x <- x$loadings
#      y <- y$loadings }
      
#  nx<- dim(x)[2]
#  ny<- dim(y)[2]
#  cross<- t(y) \%*\% x   #inner product will have dim of ny * nx
#   sumsx<- sqrt(1/diag(t(x) \%*\%x))   
 #  sumsy<- sqrt(1/diag(t(y) \%*\%y)) 

 #  result<- matrix(rep(0,nx*ny),ncol=nx)
#    result<-  sumsy * (cross * rep(sumsx, each = ny))
#  return(t(result))
#   }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
