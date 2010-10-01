\name{irt.fa}
\alias{irt.fa}

\title{Item Response Analysis by factor analysis of tetrachoric/polychoric correlations}
\description{
Although factor analysis and Item Response Theory seem to be very different models of binary data, they can provide equivalent parameter estimates of item difficulty and item discrimination.  Tetrachoric correlations of a data set of dichotomous items may be factor analysed using a minimum residual or maximum likelihood factor analysis and the result loadings transformed to item discrimination parameters.  The tau parameter from the tetrachoric correlation combined with the item factor loading may be used to estimate item difficulties. Similar analyses can be done for discrete item responses using the polychoric correlation.
}
\usage{
irt.fa(x,...)
}
\arguments{
  \item{x}{A data matrix of dichotomous or discrete items, or the result of \code{\link{tetrachoric}} or \code{\link{polychoric}}   }
  \item{...}{Additional parameters to pass to the factor analysis function}
}
\details{
The tetrachoric correlation matrix of dichotomous items may be factored using a (e.g.) minimum residual factor analyis function \code{\link{fa}} and the resulting loadings, \eqn{\lambda_i} are transformed to discriminations by
\eqn{\alpha = \frac{\lambda_i}{\sqrt{1-\lambda_i^2}} }{a = \lambda / (sqrt(1-\lambda^2)}.

The difficulty parameter, \eqn{\delta} may be found from the \eqn{\tau} parameter of the \code{\link{tetrachoric}} or \code{\link{polychoric}} function.

\eqn{\delta_i = \frac{\tau_i}{\sqrt{1-\lambda_i^2}}}{\delta = \tau / (sqrt(1-\lambda^2)}

Similar analyses may be done with discrete item responses using polychoric correlations and distinct estimates of item difficulty (location)  for each item response.

The results may be shown graphically using plot.   For plotting there are three options: type = "ICC" will plot the item characteristic respone function.  type = "IIC" will plot the item information function, and type= "test" will plot the test information function.

The normal input is just the raw data.  If, however, the correlation matrix has already been found using \code{\link{tetrachoric}} or \code{\link{polychoric}}, then that result can be processed directly.
}
\value{
 \item{coefficients }{A data frame of Item difficulty, discrimination, tau and the factor loadings}
 \item{stats }{A list of statistics for the factor analyis}
 \item{rho}{The tetrachoric/polychoric correlation matrix}

}
\references{
McDonald, Roderick P. (1999) Test theory: A unified treatment. L. Erlbaum Associates.

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{http://personality-project.org/r/book/} 
}
\author{
William Revelle
}
\note{
Still under development.  Comments welcome
}

\seealso{\code{\link{fa}}, \code{\link{sim.irt}}, \code{\link{tetrachoric}}, \code{\link{polychoric}} as well as \code{\link{plot.psych}} for plotting the IRT item curves.

}
\examples{

set.seed(17)
d9 <- sim.irt(9,1000,-2.5,2.5,mod="normal") #dichotomous items
test <- irt.fa(d9$items)
test 
op <- par(mfrow=c(3,1))
plot(test,type="ICC")
plot(test,type="IIC")
plot(test,type="test")
par(op)
set.seed(17)
items <- sim.congeneric(N=500,short=FALSE,categorical=TRUE) #500 responses to 4 discrete items
d4 <- irt.fa(items$observed)  #item response analysis of congeneric measures

op <- par(mfrow=c(2,2))
plot(d4,type="ICC")
par(op)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}
