\name{fparse}
\alias{fparse}

\title{
Parse and exten  formula input from a model and return the DV, IV, and associated terms.
}
\description{
Formula input from e.g., lm,  may be extended to include mediators,  quadratic and partial terms using a standard syntax. This is use by \code{\link{lmCor}} and \code{\link{mediate}}.
 }
\usage{
fparse(expr)
}
\arguments{
  \item{expr}{A legitimate expression in the form y ~ x1 ,  etc. (see details)}
}
\value{
\item{y}{A list of elements from the left side of the formula}
\item{x}{A list of elements from the right side of the formula}
\item{m}{A list of those elements of the formula included in ()}
\item{prod}{A list of elements separated by a * sign}
\item{ex}{A list of elements marked by I()}
}
\details{The basic formula input given as DV1 + DV2 ~ IV1 + IV2 + (IV3) + I(IV4^2) - IV5 will be parsed to return 2 DVs (1 and 2), two normal IVs (1 and 2), a mediator (IV3) a quadratic (IV4) and a variable to be partialed (IV5). See the various examples in \code{\link{lmCor}} and \code{\link{mediate}}.
}

\author{William Revelle
}



\examples{
fparse(DV  ~ IV1 + IV2 * IV2*IV3 + (IV4) + I(IV5^2) )
#somewhat more complicated
fparse(DV1 + DV2 ~ IV1 + IV2 + IV3*IV4 + I(IV5^2) + I(Iv6^2) + (IV7) + (IV8) - IV9)

}
\keyword{ utilities }
