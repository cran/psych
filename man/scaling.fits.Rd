\name{scaling.fits}
\alias{scaling.fits}

\title{ Test the adequacy of simple choice, logistic, or Thurstonian scaling. }
\description{Given a matrix of choices and a vector of scale values, how well do the scale values capture the choices?  That is, what is size of the squared residuals given the model versus the size of the squared choice values?
}
\usage{
scaling.fits(model, data, test = "logit", digits = 2, rowwise = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A vector of scale values }
  \item{data}{ A matrix or dataframe of choice frequencies }
  \item{test}{ "choice", "logistic", "normal" }
  \item{digits}{ Precision of answer }
  \item{rowwise}{Are the choices ordered by column over row (TRUE) or row over column False)}
}
\details{
}
\value{
  \item{GF }{Goodness of fit of the model}
  \item{original }{Sum of squares for original data}
  \item{resid }{Sum of squares for residuals given the data and the model}
   \item{residual }{Residual matrix}
}
\references{Revelle, W. (in preparation) Introduction to psychometric theory with applications in R, Springer.  \url{http://personality-project.org/r/book}}
\author{ William Revelle }
\note{ Mainly for demonstration purposes for a course on  psychometrics


}
\seealso{ \code{\link{thurstone}}, \code{\link{vegetables}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models}
