\name{dummy.code}
\alias{dummy.code}

\title{Create dummy coded variables}

\description{Given a variable x with n distinct values, create n new dummy coded variables coded 0/1 for presence (1) or absence (0) of each variable.  A typical application would be to create dummy coded college majors from a vector of college majors.  By default, NA values of x are returned as NA (added 10/20/17)
}
\usage{
dummy.code(x,na.rm=TRUE)
}

\arguments{
  \item{x}{A vector to be transformed into dummy codes}
  \item{na.rm}{If TRUE, return NA for all codes with NA in x}
}
\details{When coding demographic information, it is typical to create one variable with multiple categorical values (e.g., ethnicity, college major, occupation). \code{\link{dummy.code}} will convert these categories into n distinct dummy coded variables.

If using dummy coded variables as predictors, remember to use n-1 variables.}

\value{A matrix of dummy coded variables}

\author{William Revelle}


\examples{
new <- dummy.code(sat.act$education)
new.sat <- data.frame(new,sat.act)
round(cor(new.sat,use="pairwise"),2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}

