\name{unidim}
\alias{unidim}
\title{Several indices of the unidimensionality of a set of variables.}

\description{There are a variety of ways of assessing whether a set of items measures one latent trait.  \code{\link{unidim}} is just one more way.  If a one factor model holds in the data, then the factor analytic decomposition F implies that FF' should reproduce the correlations with communalities along the diagonal. In this case, the fit FF' should be identical to the correlation matrix minus the uniquenesses.  unidim is just the ratio of these two estimates.  The higher it is, the more the evidence for unidimensionality. A number of alternative statistics are estimated.

}
\usage{
unidim(x, keys.list = NULL,cor="cor",correct=.5, check.keys = TRUE)
}

\arguments{
  \item{x}{An input matrix or data frame.  If x is not a correlation matrix, then the correlations are found.}
  \item{keys.list}{If specified,  then a number of scales can be tested at once. (See \code{\link{scoreItems}} for a similar procedure.)}
\item{cor}{By default, find the Pearson correlation, other options are "spearman","kendall","tet"(for tetrachoric), "poly" (for polychoric), or  "mixed"}
\item{correct}{If using "tetrachoric" or "polychoric" correlations, should we correct empty cells for continuity, and if so, by how much.  (See \code{\link{tetrachoric}} for a discussion of this correction)}
  \item{check.keys}{If TRUE, then items will be keyed based upon their loadings on the first factor.  Automatically done if key.list is NULL.}

}
\details{
This is set of exploratory indices that are still under development.  A number of test cases suggest that u provides high values when the data are in fact unidimensional, low values when they are not.

The logic is deceptively simple:  Unidimensionality implies that a one factor model of the data fits the covariances of the data.  If this is the case, then factor model implies R = FF' + U2 will have residuals of 0.  Similarly, this also implies that the observed correlations will equal the model.  Thus, the sum of the observed correlations (with the diagonal replaced by the communalities) should match the factor model.  Compare these two models:  R - U2  versus FF'.  This is the unidim.A estimate. 

Also reported, and probably better, is the fit of the one factor model to the correlations. This is merely the sumsq of the residual correlations/sumsq of the original correlations.  When the factor model is perfect, this will be 1. 

This works well, but when some of the loadings are very small, even though 1 factor is correct, it is probably not a good idea to think of the items as forming a unidimensional scale.  

An alternative model (the av.r.fit statistic) considers the residuals found by subtracting the average correlation from the observed correlations.  This will achieve a maximum if the item covariances are all identical (a tau equivalent model).

The product of fa.fit and av.r.fit is the measure of unidimensionality, u.  
 
 All of the results are reported in the uni list.
}
\value{ 
\item{u}{The estimate of unidimensionality which is just the product of}
\item{av.r.fit}{The fit of the average r to the correlation matrix}
\item{fa.fit}{The off diagonal fit from  \code{\link{fa}}  }
\item{alpha}{Standardized alpha of the keyed items (after appropriate reversals)}
\item{av.r}{The average interitem correlation of the keyed items.}
\item{median.r}{The median value of the iteritem correlations of the keyed items.}

\item{Unidim.A}{The unidimensional criterion when items are keyed in positive direction.}

\item{Unidim }{The raw value of the unidimensional criterion}
\item{raw.model}{The ratio of the FF' model to the sum(R)}
\item{adj.model}{The ratio of the FF' model to the sum(R) when items are flipped.}
\item{Total}{The ratio of the sum(R - uniqueness)/sum(R)}
\item{Total.A}{Same ratio with flipped items}


}

\author{William Revelle}

\note{A perhaps interesting idea but still an exploratory statistic. Treat with appropriate caution. 
}


\seealso{ \code{\link{fa}} for factor analysis, \code{\link{omega}} for reliability.

}
\examples{
#test the unidimensionality of the five factors of the bfi data set.


unidim(psychTools::bfi,psychTools::bfi.keys) 
unidim(psychTools::ability,psychTools::ability.keys)
#Try a known 3 factor structure
x <- sim.minor(nfact=3,bipolar=FALSE)   #this makes all the items positive
unidim(x$model) 
keys.list <- list(first =c(1:4),second = 5:8,third=9:12,all=1:12)
unidim(x$model,keys.list)

x <- sim.minor(nfact=3)
unidim(x$model,keys.list)   #we flip the negative items 

#what about a hierarchical model?
H <- sim.hierarchical()  # by default, a nice hierarchical model
H.keys <- list(First = paste0("V",1:3),Second=paste0("V",4:6),Third=paste0("V",7:9),
  All = paste0("V",1:9))
unidim(H,H.keys)



}


\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line