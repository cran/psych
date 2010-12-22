\name{irt.fa}
\alias{irt.fa}
\alias{irt.select}

\title{Item Response Analysis by factor analysis of tetrachoric/polychoric correlations}
\description{
Although factor analysis and Item Response Theory seem to be very different models of binary data, they can provide equivalent parameter estimates of item difficulty and item discrimination.  Tetrachoric correlations of a data set of dichotomous items may be factor analysed using a minimum residual or maximum likelihood factor analysis and the result loadings transformed to item discrimination parameters.  The tau parameter from the tetrachoric correlation combined with the item factor loading may be used to estimate item difficulties. Similar analyses can be done for discrete item responses using the polychoric correlation.
}
\usage{
irt.fa(x,...)
irt.select(x,y)
}
\arguments{
  \item{x}{A data matrix of dichotomous or discrete items, or the result of \code{\link{tetrachoric}} or \code{\link{polychoric}}   }
  \item{y}{the subset of variables to pick from the rho and tau output of a previous irt.fa analysis to allow for further analysis.}
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

The output is best seen in terms of graphic displays.  Plot the output from irt.fa to see item and test information functions.  

The select.irt function is a helper function to allow for selecting a subset of a prior analysis for further analysis. First run irt.fa, then select a subset of variables to be analyzed in a subsequent irt.fa analysis.
}
\value{
 \item{irt }{A list of Item location (difficulty) and discrimination}
 \item{fa }{A list of statistics for the factor analyis}
\item{rho}{The tetrachoric/polychoric correlation matrix}
\item{tau}{The tetrachoric/polychoric cut points}

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


#using the iq data set for an example of real items
#first need to convert the responses to tf
data(iqitems)
iq.keys <- c(4,4,3,1,4,3,2,3,1,4,1,3,4,3)
iq.tf <- score.multiple.choice(iq.keys,iqitems,score=FALSE)  #just the responses
iq.irt <- irt.fa(iq.tf)
plot(iq.irt)
#select a subset of these variables
small.iq.irt <- irt.select(iq.irt,c(1,5,9,10,11,13))
small.irt <- irt.fa(small.iq.irt)
plot(small.irt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}
