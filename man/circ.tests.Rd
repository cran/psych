\name{circ.tests}
\alias{circ.tests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Apply four tests of circumplex versus simple structure }
\description{Rotations of factor analysis and principal components analysis solutions typically try to represent correlation matrices as simple structured.  An alternative structure, appealing to some, is a circumplex structure where the variables are uniformly spaced on the perimeter of a circle in a two dimensional space.  Generating these data is straightforward, and is useful for exploring alternative solutions to affect and personality structure.
}
\usage{
circ.tests(loads, loading = TRUE, sorting = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loads}{ A matrix of loadings \code{loads} here }
  \item{loading}{ Are these loadings or a correlation matrix \code{loading} }
  \item{sorting}{ Should the variables be sorted \code{sorting} }
}
\details{``A common model for representing psychological data is simple structure (Thurstone, 1947). According to one common interpretation, data are simple structured when items or scales have non-zero factor loadings on one and only one factor (Revelle & Rocklin, 1979). Despite the commonplace application of simple structure, some psychological models are defined by a lack of simple structure. Circumplexes (Guttman, 1954) are one kind of model in which simple structure is lacking.

``A number of elementary requirements can be teased out of the idea of circumplex structure. First, circumplex structure implies minimally that variables are interrelated; random noise does not a circumplex make. Second, circumplex structure implies that the domain in question is optimally represented by two and only two dimensions. Third, circumplex structure implies that variables do not group or clump along the two axes, as in simple structure, but rather that there are always interstitial variables between any orthogonal pair of axes (Saucier, 1992). In the ideal case, this quality will be reflected in equal spacing of variables along the circumference of the circle (Gurtman, 1994; Wiggins, Steiger, & Gaelick, 1981). Fourth, circumplex structure implies that variables have a constant radius from the center of the circle, which implies that all variables have equal communality on the two circumplex dimensions (Fisher, 1997; Gurtman, 1994). Fifth, circumplex structure implies that all rotations are equally good representations of the domain (Conte & Plutchik, 1981; Larsen & Diener, 1992). (Acton and Revelle, 2004)

Acton and Revelle reviewed the effectiveness of 10 tests of circumplex structure and found that four did a particularly good job of discriminating circumplex structure from simple structure, or circumplexes from ellipsoidal structures. Unfortunately, their work was done in Pascal and is not easily available. Here we release R code to do the four most useful tests:

	1	The Gap test of equal spacing
	
	2	Fisher's test of equality of axes 
	
	3	A test of indifference to Rotation
	
	4	A test of equal Variance of squared factor loadings across arbitrary rotations.
	
To interpret the values of these various tests, it is useful to compare the particular solution to simulated solutions representing pure cases of circumplex and simple structure.  See the example output from \code{\link{circ.simulation}} and compare these plots with the results of the circ.test.
}
\value{A list of four items is returned.  These are the gap, fisher, rotation and variance test results.
  \item{gaps}{gap.test}
  \item{fisher}{fisher.test}
  \item{RT}{rotation.test}
  \item{VT}{variance.test}
 
}
\references{ Acton, G. S. and Revelle, W. (2004) Evaluation of Ten Psychometric Criteria for Circumplex Structure.  Methods of Psychological Research Online, Vol. 9, No. 1 \url{http://www.dgps.de/fachgruppen/methoden/mpr-online/issue22/mpr110_10.pdf} }
\author{ William Revelle}
\note{ Of the 10 criterion discussed in Acton and Revelle (2004), these tests operationalize the four most useful.

}
\seealso{  \code{\link{circ.simulation}}, \code{\link{circ.sim}} }
\examples{
circ.data <- circ.sim(24,500)
circ.fa <- factor.pa(circ.data,2)
#plot(circ.fa$loadings)
ct <- circ.tests(circ.fa)
#compare with non-circumplex data
simp.data <- item.sim(24,500)
simp.fa <- factor.pa(simp.data,2)
#plot(simp.fa$loadings)
st <- circ.tests(simp.fa)
print(rbind(ct,st),digits=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}

