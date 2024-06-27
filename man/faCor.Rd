\name{faCor}
\alias{faCor}
\title{Correlations between two factor analysis solutions}
\description{
Given two factor analysis or pca solutions to a data matrix or correlation, what are the similarities between the two solutions. This may be found by factor correlations as well as factor congruences.  Factor correlations are found by the matrix product of the factor weights and the correlation matrix and are estimates of what the factor score correlations would be.  Factor congruence (aka Tucker or Burt coefficient) is the cosine of the vectors of factor loadings.
}
\usage{
faCor(r, nfactors = c(1, 1), fm = c("minres", "minres"), rotate =
 c("oblimin", "oblimin"), scores = c("tenBerge", "tenBerge"), adjust=c(TRUE,TRUE),
     use = "pairwise", cor = "cor", weight = NULL, correct = 0.5,Target=list(NULL,NULL))
}
\arguments{
  \item{r}{A correlation matrix or a data matrix suitable for factoring}
  \item{nfactors}{Number of factors in each solution to extract}

  \item{fm}{Factor method.  The default is 'minres' factoring. To compare with pca solutions, can also be (fm ="pca")
}
  \item{rotate}{What type of rotations to apply.  The default for factors is oblimin, for pca is varimax.}
 \item{scores}{What factor scoring algorithm should be used. Defaults to tenBerge for both cases.}
 \item{adjust}{Should the factor intercorrelations be corrected by the lack of factor deteriminancy (i.e. divide by the square root of the factor R2)}
 
  \item{use}{How to treat missing data.  Use='pairwise" finds pairwise complete correlations.
}
  \item{cor}{What kind of correlation to find.  The default is Pearson.
}
  \item{weight}{Should cases be weighted?  Default, no.}
  \item{correct}{If finding tetrachoric or polychoric correlations, what correction should be applied to empty cells (defaults to .5)}
\item{Target}{If doing target rotations (e.g., TargetQ or TargetT), then the Target must be specified.  If TargetT, this may be a matrix, if TargetQ, this must be a list.  (Strange property of GPARotation.)}
 
}

\details{The factor correlations are found using the approach discussed by Gorsuch (1983) and uses the weights matrices found by \eqn{w=S R^{-1}} and \eqn{r = w' R w} where S is the structure matrix and is   \eqn{s= F \Phi}.  The resulting correlations may be  adjusted for the factor score variances (the diagonal of r) (the default). 

For factor loading vectors of F1 and F2 the measure of factor congruence, phi, is 
\deqn{
\phi = \frac{\sum F_1 F_2}{\sqrt{\sum(F_1^2)\sum(F_2^2)}} 
.}{{phi = sum(F1 F2)/sqrt(sum(F1^2) sum(F2^2))} } and is also found in \code{\link{factor.congruence}}.

For comparisons of factor solutions from 1 to n, use \code{\link{bassAckward}}. This function just compares two solutions from the same correlation/data matrix.  \code{\link{factor.congruence}} can be used to compare any two sets of factor loadings. 

Note that alternative ways of finding weights (e.g., regression, Bartlett, tenBerge) will produce somewhat different results.  tenBerge produces weights that maintain the factor correlations in the factor scores.
}
\value{
\item{Call}{The function call}
\item{r}{The factor intercorrelations}
\item{congruence}{The Burt/Tucker coefficient of congruence}
\item{f1}{The first factor analysis}
\item{f2}{The second factor analysis}
}
\references{

Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates.

Burt, Cyril (1948) The factorial study of temperamental traits. British Journal of Statistical Psychology, 1(3) 178-203.

 Horn, John L. (1973) On extension analysis and its relation to correlations between variables and factor scores.  Multivariate Behavioral Research, 8, (4), 477-489.

Lorenzo-Seva, U. and ten Berge, J. M. F. (2006). Tucker's congruence coefficient as a meaningful index of factor similarity. Methodology: European Journal of Research Methods for the Behavioral and Social Sciences, 2(2):57-64.

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 
}
\author{
William Revelle
}
\note{Useful for comparing factor solutions from the same data.  Will not work for different data sets
}


\seealso{ \code{\link{fa}}, \code{\link{pca}},  \code{\link{omega}} and \code{\link{iclust}}, and \code{\{link{bassAckward}} for alternative hierarchical solutions.  \code{\link{fa.extend}} and \code{\link{fa.extension}} for other uses of factor - item correlations.
}
\examples{
faCor(Thurstone,nfactors=c(2,3)) #compare two solutions to the Thurstone problem
faCor(bfi[1:25],nfactors=c(5,5),fm=c("minres","pca")) #compare pca and fa solutions
#compare two levels of factor extraction, and then find the correlations of the scores
faCor(bfi[1:25],nfactors=c(3,5)) #based upon linear algebra 
f3 <- fa(bfi[1:25],3,scores="tenBerge")
f5 <- fa(bfi[1:25],5 ,scores="tenBerge")
cor2(f3$scores,f5$scores) #the correlation between the factor score estimates
  
}
\keyword{multivariate }
\keyword{ models}