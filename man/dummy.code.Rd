\name{dummy.code}
\alias{dummy.code}

\title{Create dummy coded variables}

\description{Given a variable x with n distinct values, create n new dummy coded variables coded 0/1 for presence (1) or absence (0) of each variable.  A typical application would be to create dummy coded college majors from a vector of college majors. Can also combine categories by group.  By default, NA values of x are returned as NA (added 10/20/17)
}
\usage{
dummy.code(x,group=NULL,na.rm=TRUE)
}

\arguments{
  \item{x}{A vector to be transformed into dummy codes}
  \item{group}{A vector of categories to be coded as 1, all others coded as 0.}
  \item{na.rm}{If TRUE, return NA for all codes with NA in x}
}
\details{When coding demographic information, it is typical to create one variable with multiple categorical values (e.g., ethnicity, college major, occupation). \code{\link{dummy.code}} will convert these categories into n distinct dummy coded variables.

If using dummy coded variables as predictors, remember to use n-1 variables.

If group is specified, then all values of x that are in group are given the value of 1, otherwise, 0. (Useful for combining a range of science majors into STEM or not. The example forms a dummy code of any smoking at all.)
}

\value{A matrix of dummy coded variables}

\author{William Revelle}


\examples{
new <- dummy.code(sat.act$education)
new.sat <- data.frame(new,sat.act)
round(cor(new.sat,use="pairwise"),2)
dum.smoke <- dummy.code(spi$smoke,group=2:9)
table(dum.smoke,spi$smoke)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}

