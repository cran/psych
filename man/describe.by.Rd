\name{describe.by}
\alias{describe.by}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Basic summary statistics by group}
\description{Report basic summary statistics by a grouping variable.  Useful if the grouping variable is some experimental variable and data are to be aggregated for plotting.  Just a wrapper for by and \code{\link{describe}}.
}
\usage{
describe.by(x, group,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data.frame or matrix  }
  \item{group}{a grouping variable}
  \item{...}{parameters to be passed to describe}
}
\details{
}
\value{
 A data.frame of the relevant statistics broken down by group: \cr
   item name \cr
   item number \cr
   number of valid cases\cr
   mean\cr
   standard deviation\cr
   median\cr
   mad: median absolute deviation (from the median) \cr
   minimum\cr
   maximum\cr
   skew\cr
   standard error\cr
}

\author{ William Revelle}

\seealso{  \code{\link{describe}}}
\examples{
set.seed(42)  #to get the same values each time
 x.df <- data.frame(group=sample(2,20,replace=TRUE), matrix(rnorm(100),ncol=5))
 x <- describe.by(x.df,x.df$group)
 x  #shows all the results
 x[1]  #shows just the first group
 x <- matrix(sample(4,200,replace=TRUE),ncol=5)
 y <- describe.by(x,x[,1])
y


}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ univar }% __ONLY ONE__ keyword per line
