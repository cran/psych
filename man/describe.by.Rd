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
#
#group: 1
#      var n  mean   sd median  mad   min  max range   se
#group   1 7  1.00 0.00   1.00 0.00  1.00 1.00  0.00 0.00
#X1      2 7 -0.36 1.66  -0.31 2.16 -2.66 1.90  4.55 0.63
#X2      3 7 -0.10 1.04  -0.43 0.95 -1.37 1.44  2.81 0.39
#...
#----------------------------------------------------------------------------- 
#group: 2
#      var  n  mean   sd median  mad   min  max range   se
#group   1 13  2.00 0.00   2.00 0.00  2.00 2.00  0.00 0.00
#X1      2 13 -0.07 1.26  -0.26 0.57 -2.44 2.29  4.73 0.35
#...
#>  x[1]  #shows just the first group
#$`1`
#      var n  mean   sd median  mad   min  max range   se
#group   1 7  1.00 0.00   1.00 0.00  1.00 1.00  0.00 0.00
#X1      2 7 -0.36 1.66  -0.31 2.16 -2.66 1.90  4.55 0.63
#X2      3 7 -0.10 1.04  -0.43 0.95 -1.37 1.44  2.81 0.39
#X3      4 7  0.56 0.90   0.92 0.97 -0.73 1.58  2.30 0.34
#X4      5 7  0.29 0.88   0.46 0.87 -1.19 1.51  2.71 0.33
#X5      6 7  0.23 1.12   0.19 1.27 -1.46 1.85  3.31 0.42


## The function is currently defined as
function (x,group,...) {               #data are x, grouping variable is group
answer <- by(x,group,describe,...)
return(answer)}
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
