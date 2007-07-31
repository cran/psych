\name{test.psych}
\alias{test.psych}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Testing of functions in the psych package }
\description{Test to make sure the psych functions run on basic test data sets
}
\usage{
test.psych(first=1,last=5,short=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{first}{first=1: start with dataset first}
  \item{last}{last=5: test for datasets until last}
  \item{short}{short=TRUE - don't return any analyses}
}
\details{When modifying the psych package, it is useful to make sure that adding some code does not break something else.  The test.psych function tests the major functions on various standard data sets.  It  also shows off a number of the capabilities of the psych package.
}
\value{
  \item{out }{if short=FALSE, then list of the output from all functions tested}

}

\author{ William Revelle}
\note{ Not for general consumption -- used to make sure functions throw error messages or correct for weird conditions.

The datasets tested are part of the standard R data sets and represent some of the basic problems encountered.}

\section{Warning }{Warning messages will be thrown by fa.parallel and sometimes by factor.pa for random datasets.} 
}

\examples{
#test <- psych.test()


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

