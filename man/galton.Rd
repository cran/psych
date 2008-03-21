\name{galton}
\alias{galton}
\docType{data}
\title{Galton's Mid parent child height data}
\description{Two of the earliest examples of the correlation coefficient were Francis Galton's data sets on the relationship between mid parent and child height and the similarity of parent generation peas with child peas.  This is the data set for the Galton height.
}
\usage{data(galton)}
\format{
  A data frame with 928 observations on the following 2 variables.
  \describe{
    \item{\code{parent}}{Mid Parent heights (in inches) }
    \item{\code{child}}{Child Height}
  }
}
\details{Female heights were adjusted by 1.08 to compensate for sex differences. (This was done in the original data set)
}
\source{This is just the galton data set from UsingR, slightly rearranged.
}
\references{Stigler, S. M. (1999). Statistics on the Table: The History of Statistical 
Concepts and Methods. Harvard University Press.       
  Galton, F. (1869). Hereditary Genius: An Inquiry into its Laws and Con- 
sequences. London: Macmillan. 
 
  Wachsmuth, A.W., Wilkinson L., Dallal G.E. (2003). Galton's bend: A 
previously undiscovered nonlinearity in Galton's family stature regression 
data. The American Statistician, 57, 190-1922.
}
\examples{
data(galton)
describe(galton)
pairs.panels(galton,lm=TRUE)
}
\keyword{datasets}
