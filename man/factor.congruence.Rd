\name{factor.congruence}
\alias{factor.congruence}
\alias{fa.congruence}
\title{Coefficient of factor congruence }
\description{Given two sets of factor loadings, report their degree of congruence (vector cosine). Although reported by Burt (1948), this is frequently known as the Tucker index of factor congruence.
}
\usage{
factor.congruence(x, y=NULL,digits=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix of factor loadings or a list of matrices of factor loadings}
  \item{y}{ A second matrix of factor loadings (if x is a list, then y may be empty)}
  \item{digits}{Round off to digits}

}
\details{Find the coefficient of factor congruence between two sets of factor loadings. 

Factor congruences are the cosines of pairs of vectors defined by the loadings matrix and based at the origin.  Thus, for loadings that differ only by a scaler (e.g. the size of the eigen value), the factor congruences will be 1.

For factor loading vectors of F1 and F2 the measure of factor congruence, phi, is 
\deqn{
\phi = \frac{\sum F_1 F_2}{\sqrt{\sum(F_1^2)\sum(F_2^2)}} 
.}{phi = sum(F1 F2)/sqrt(sum(F1^2) sum(F2^2)) }

It is an interesting exercise to compare factor congruences with the correlations of factor loadings.  Factor congruences are based upon the raw cross products, while correlations are based upon centered cross products. That is,  correlations of factor loadings are cosines of the vectors based at the mean loading for each factor.   

\deqn{
\phi = \frac{\sum (F_1-a) (F_2 - b)}{\sqrt{\sum((F_1-a)^2)\sum((F_2-b)^2)}} 
.}{phi = sum((F1-a)(F2-b))/sqrt(sum((F1-a)^2) sum((F2-b)^2)) }.

For congruence coefficients, a = b= 0.  For correlations a=mean F1, b= mean F2.

Input may either be matrices or factor analysis or principal components analyis output (which includes a loadings object), or a mixture of the two.

To compare more than two solutions, x may be a list of matrices, all of which will be compared.


}
\value{A matrix of factor congruences.
}

\references{
Burt, Cyril (1948) The factorial study of temperamental traits. British Journal of Statistical Psychology, 1(3) 178-203.

Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates.
\cr

Revelle, W. (In preparation) An Introduction to Psychometric Theory with applications in R (\url{http://personality-project.org/r/book/})
}
\author{
 \email{revelle@northwestern.edu} \cr
 \url{http://personality-project.org/revelle.html}}

\seealso{ \code{\link{principal}}, \code{\link{fa}}}
 
\examples{
#factor congruence of factors and components, both rotated
#fa <- fa(Harman74.cor$cov,4)
#pc <- principal(Harman74.cor$cov,4)
#factor.congruence(fa,pc)
 #    RC1  RC3  RC2  RC4
#MR1 0.98 0.41 0.28 0.32
#MR3 0.35 0.96 0.41 0.31
#MR2 0.23 0.16 0.95 0.28
#MR4 0.28 0.38 0.36 0.98



#factor congruence without rotation
#fa <- fa(Harman74.cor$cov,4,rotate="none")
#pc <- principal(Harman74.cor$cov,4,rotate="none")
#factor.congruence(fa,pc)   #just show the beween method congruences
#     PC1   PC2   PC3   PC4
#MR1 1.00 -0.04 -0.06 -0.01
#MR2 0.15  0.97 -0.01 -0.15
#MR3 0.31  0.05  0.94  0.11
#MR4 0.07  0.21 -0.12  0.96

#factor.congruence(list(fa,pc))  #this shows the within method congruence as well

 #     MR1   MR2  MR3   MR4  PC1   PC2   PC3   PC4
#MR1  1.00  0.11 0.25  0.06 1.00 -0.04 -0.06 -0.01
#MR2  0.11  1.00 0.06  0.07 0.15  0.97 -0.01 -0.15
#MR3  0.25  0.06 1.00  0.01 0.31  0.05  0.94  0.11
#MR4  0.06  0.07 0.01  1.00 0.07  0.21 -0.12  0.96
#PC1  1.00  0.15 0.31  0.07 1.00  0.00  0.00  0.00
#PC2 -0.04  0.97 0.05  0.21 0.00  1.00  0.00  0.00
#PC3 -0.06 -0.01 0.94 -0.12 0.00  0.00  1.00  0.00
#PC4 -0.01 -0.15 0.11  0.96 0.00  0.00  0.00  1.00

#pa <- fa(Harman74.cor$cov,4,fm="pa")
# factor.congruence(fa,pa)
#         PA1  PA3  PA2  PA4
#Factor1 1.00 0.61 0.46 0.55
#Factor2 0.61 1.00 0.50 0.60
#Factor3 0.46 0.50 1.00 0.57
#Factor4 0.56 0.62 0.58 1.00


#compare with 
#round(cor(fa$loading,pc$loading),2)
#      RC1   RC3   RC2   RC4
#MR1  0.99 -0.18 -0.33 -0.34
#MR3 -0.33  0.96 -0.16 -0.43
#MR2 -0.29 -0.46  0.98 -0.21
#MR4 -0.44 -0.30 -0.22  0.98

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
