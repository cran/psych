\name{polar}
\alias{polar}
\title{Convert Cartesian factor loadings into polar coordinates }
\description{Factor and cluster analysis output typically presents item by factor correlations (loadings).  Tables of factor loadings are frequently sorted by the size of loadings.  This style of presentation tends to make it difficult to notice the pattern of loadings on other, secondary, dimensions.  By converting to polar coordinates, it is easier to see the pattern of the secondary loadings. 
}
\usage{
polar(f, sort = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{A matrix of loadings or the output from a factor or cluster analysis program}
  \item{sort}{ sort=TRUE: sort items by the angle of the items on the first pair of factors.}
}
\details{Although many uses of factor analysis/cluster analysis assume a simple structure where items have one and only one large loading, some domains such as personality or affect items have a more complex structure and some items have high loadings on two factors.  (These items are said to have complexity 2, see \code{\link{VSS}}).  By expressing the factor loadings in polar coordinates, this structure is more readily perceived.

For each pair of factors, item loadings are converted to an angle with the first factor, and a vector length corresponding to the amount of variance in the item shared with the two factors.  

For a two dimensional structure, this will lead to a column of angles and a column of vector lengths.  For n factors, this leads to n* (n-1)/2 columns of angles and an equivalent number of vector lengths.

}
\value{

  \item{polar }{A data frame of polar coordinates }
 
}
\references{Rafaeli, E. & Revelle, W. (2006). A premature consensus: Are happiness and sadness truly opposite affects? Motivation and Emotion. \\

Hofstee, W. K. B., de Raad, B., & Goldberg, L. R. (1992). Integration of the big five and circumplex approaches to trait structure. Journal of Personality and Social Psychology, 63, 146-163.}
\author{William Revelle}

\seealso{ \code{\link{ICLUST}}, \code{\link{cluster.plot}}, \code{\link{circ.tests}},  \code{\link{factor.pa}} }
\examples{

circ.data <- circ.sim(24,500)
circ.fa <- fa(circ.data,2)
circ.polar <- round(polar(circ.fa),2)
circ.polar
#compare to the graphic
cluster.plot(circ.fa)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

