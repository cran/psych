\name{alpha.scale}
\alias{alpha.scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cronbach alpha for a scale }
\description{
  Find Cronbach's coefficient alpha given a scale and a data.frame of the items in the scale.
  For X, a total score composed of items in the data.frame Y, find Cronbach's alpha.
}
\usage{
alpha.scale(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A scale made by summing items }
  \item{y}{ A data frame of items }
}

\value{ Coefficient alpha }
\details{Alpha is one of several estimates of the internal consistency reliability of a test.
Perhaps because it is so easy to calculate, it is without doubt the most frequently reported measure of internal consistency reliability. Alpha is the mean of all possible spit half reliabilities (corrected for test length).  For a unifactorial test, it is a reasonable estimate of the first factor saturation, although if the test has any microstructure (i.e., if it is ``lumpy") coefficients beta (see \code{\link{ICLUST}} and \code{link{omega}} are more appropriate estimates of the general factor saturation.

Alpha is a positive function of the number of items in a test as well as the average intercorrelation of the items in the test.  When calculated from the item variances and total test variance, as is done here, alpha is sensitive to differences in the item variances.  Alternative functions \code{\link{score.items}} and  \code{\link{cluster.cor}} will also score multiple scales and report more useful statistics. ``Standardized" alpha is calculated from the inter-item correlations and will differ from raw alpha. Standardized alpha can be found by using \code{\link{cluster.cor}}. }


\references{ \url{http://personality-project.org/revelle/syllabi/405.syllabus.html}}
\author{

Maintainer: William Revelle \email{revelle@northwestern.edu}
}


\seealso{\code{\link{score.items}}, \code{\link{cluster.cor}}, \code{\link{ICLUST}}, \code{\link{omega}}, \code{\link{describe}},\code{\link{pairs.panels}}, alpha in psychometrics }
\examples{

y <- attitude     #from the datasets package
x <- rowSums(y)   #find the sum of the seven attitudes
alpha.scale(x,y)
#[1] 0.8431428
#compare with standardized alpha:
st.alpha <- cluster.cor(rep(1,7),cor(attitude),digits=4)
st.alpha 
#compare with score.items
si <- score.items(rep(1,7), attitude,digits=3)
si$alpha


}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line
