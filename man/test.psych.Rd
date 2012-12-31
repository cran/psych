\name{test.psych}
\alias{test.psych}
\title{ Testing of functions in the psych package }
\description{Test to make sure the psych functions run on basic test data sets
}
\usage{
test.psych(first=1,last=5,short=TRUE,all=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{first}{first=1: start with dataset first}
  \item{last}{last=5: test for datasets until last}
  \item{short}{short=TRUE - don't return any analyses}
  \item{all}{To get around a failure on certain Solaris 32 bit systems, all=FALSE is the default}
}
\details{When modifying the psych package, it is useful to make sure that adding some code does not break something else.  The test.psych function tests the major functions on various standard data sets.  It  also shows off a number of the capabilities of the psych package.

Uses 5 standard data sets: \cr
USArrests               Violent Crime Rates by US State  (4 variables) \cr
attitude               The Chatterjee-Price Attitude Data    \cr
Harman23.cor\$cov       Harman Example 2.3 8 physical measurements  \cr
Harman74.cor\$cov        Harman Example 7.4  24 mental measurements \cr
ability.cov\$cov        8 Ability and Intelligence Tests  \cr
 
}
\value{
  \item{out }{if short=FALSE, then list of the output from all functions tested}

}

\author{ William Revelle}
\note{ Although test.psych may be used as a quick demo of the various functions in the psych packge, in general, it is better to try the specific functions themselves.  The main purpose of test.psych is to make sure functions throw error messages or correct for weird conditions.

The datasets tested are part of the standard R data sets and represent some of the basic problems encountered.

When version 1.1.10 was released, it caused errors when compiling and testing on some Solaris 32 bit systems.  The all option was added to avoid this problem (since I can't replicate the problem on Macs or PCs).  all=TRUE adds one more test, for a non-positive definite matrix.}

\section{Warning }{Warning messages will be thrown by fa.parallel and sometimes by factor.pa for random datasets.} 


\examples{
#test <- test.psych()
#not run
#test.psych(all=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

