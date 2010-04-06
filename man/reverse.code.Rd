\name{reverse.code}
\alias{reverse.code}

\title{Reverse the coding of selected items prior to scale analysis}

\description{Some IRT functions require all items to be coded in the same direction.  Some data sets have items that need to be reverse coded (e.g., 6 -> 1, 1 -> 6).  reverse.code will flip items based upon a keys vector of 1s and -1s.  Reversed items are subtracted from the item max + item min.  These may be specified or may be calculated.
}
\usage{
reverse.code(keys, items, mini = NULL, maxi = NULL)}

\arguments{
  \item{keys}{A vector of 1s and -1s.  -1 implies reverse the item}
  \item{items}{A data set of items}
  \item{mini}{if NULL, the empirical minimum for each item.  Otherwise, a vector of minima}
  \item{maxi}{f NULL, the empirical maximum for each item.  Otherwise, a vector of maxima}
}
\details{Not a very complicated function, but useful in the case that items need to be reversed prior to using IRT functions from the ltm or eRM packages. Most psych functions do not require reversing prior to analysis, but will do so within the function.
}
\value{
The corrected items.
}

\examples{
original <- matrix(sample(6,50,replace=TRUE),10,5)
keys <- c(1,1,-1,-1,1)  #reverse the 3rd and 4th items
new <- reverse.code(keys,original,mini=rep(1,5),maxi=rep(6,5))
original[1:3,]
new[1:3,]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
