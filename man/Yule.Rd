\name{Yule}
\alias{Yule}
\alias{Yule.inv}
\alias{Yule2phi}
\alias{Yule2poly}
\title{From a two by two table, find the Yule coefficients of association, convert to phi, or polychoric, recreate table the table to create the Yule coefficient.}
\description{One of the many measures of association is the Yule coefficient.  Given a two x two table of counts \cr
\tabular{lll}{
\tab a \tab b \cr
\tab c \tab d  \cr
}
Yule Q is (ad - bc)/(ad+bc). \cr
Conceptually, this is the number of pairs in agreement (ad) - the number in disagreement (bc) over the total number of paired observations.  
\cr
ad/bc is the odds ratio and Q = (OR-1)/(OR+1) 
\cr
Yule's coefficient of colligation is Y = (sqrt(OR) - 1)/(sqrt(OR)+1)
Yule.inv finds the cell entries for a particular Q and the marginals (a+b,c+d,a+c, b+d).  This is useful for converting old tables of correlations into more conventional \code{\link{phi}} or polychoric correlations.
\cr
Yule2phi and Yule2poly convert the Yule Q with set marginals to the correponding phi or polychoric correlation.
}
\usage{
Yule(x,Y=FALSE)
Yule.inv(Q,m)
Yule2phi(Q,m)
Yule2poly(Q,m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of four elements or a two by two matrix }
  \item{Y}{Y=TRUE return Yule's Y coefficient of colligation}
  \item{Q}{The Yule coefficient}
  \item{m}{A two x two matrix of marginals or a four element vector of marginals}
}
\details{Yule developed two measures of association for two by two tables.  Both are functions of the odds ratio 
}
\value{
  \item{Q}{The Yule Q coefficient}
  \item{R}{A two by two matrix of counts}
}
\references{Yule, G. Uday (1912) On the methods of measuring association between two attributes. Journal of the Royal Statistical Society, LXXV, 579-652}

\author{ William Revelle }
\note{Currently done by using the optimize function, but presumably could be redone by solving a quadratic equation.
}
\seealso{ See Also as \code{\link{phi}}, \code{\link{Yule2poly.matrix}}, \code{\link{Yule2phi.matrix}} }
\examples{
Nach <- matrix(c(40,10,20,50),ncol=2,byrow=TRUE)
Yule(Nach)
Yule.inv(.81818,c(50,70,60,60))
Yule2phi(.81818,c(50,70,60,60))
if(require(polycor)) Yule2poly(.81818,c(50,70,60,60))
phi(Nach)  #much less
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{models }% __ONLY ONE__ keyword per line