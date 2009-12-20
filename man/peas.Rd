\name{peas}
\alias{peas}
\docType{data}
\title{Galton`s Peas}
\description{Francis Galton introduced the correlation coefficient with an analysis of the similarities of the parent and child generation of 700 sweet peas.
}
\usage{data(peas)}
\format{
  A data frame with 700 observations on the following 2 variables.
  \describe{
    \item{\code{parent}}{The mean diameter of the mother pea for 700 peas}
    \item{\code{child}}{The mean diameter of the daughter pea for 700 sweet peas}
  }
}
\details{Galton's introduction of the correlation coefficient was perhaps the most important contribution to the study of individual differences.  This data set allows a graphical analysis of the data set.  There are two different graphic examples.  One shows the regression lines for both relationships, the other finds the correlation as well.
}
\source{Stanton, Jeffrey M. (2001) Galton, Pearson, and the Peas: A brief history of linear regression for statistics intstructors, Journal of Statistics Education, 9. (retrieved from the web from http://www.amstat.org/publications/jse/v9n3/stanton.html) reproduces the table from Galton, 1894, Table 2. 

The data were generated from this table.
}
\references{Galton, Francis (1877) Typical laws of heredity. paper presented to the weekly evening meeting of the Royal Institution, London. Volume VIII (66) is the first reference to this data set.  The data appear in 

Galton, Francis (1894) Natural Inheritance (5th Edition), New York: MacMillan). 
}

\seealso{The other Galton data sets:   \code{\link{heights}}, \code{\link{galton}},\code{\link{cubits}}}
\examples{
data(peas)
pairs.panels(peas,lm=TRUE,xlim=c(14,22),ylim=c(14,22),main="Galton's Peas")
describe(peas)
pairs.panels(peas,main="Galton's Peas")
}
\keyword{datasets}
