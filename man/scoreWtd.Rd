\name{scoreWtd}
\alias{scoreWtd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Score items using regression or correlation based weights}
\description{Item weights from \code{\link{bestScales}} or \code{\link{setCor}} are used to find weighted scale scores. In contrast to the unit weights used in \code{\link{scoreItems}}, \code{\link{scoreWtd}} will multiply the data by a set of weights to find scale scores.  These weight may come from a regression (e.g., \code{\link{lm}} or \code{\link{setCor}}) or may be the zero order correlation weights from \code{\link{bestScales}}.
}
\usage{
scoreWtd(weights, items, std = TRUE, sums = FALSE, impute = "none")
}

\arguments{
\item{weights}{This is just a matrix of weights to use for each item for each scale.}
\item{items}{ Matrix or dataframe of raw item scores}
\item{std}{if TRUE, then find weighted standard scores else just use raw data}
\item{sums}{By default, find the average item score. If sums = TRUE, then find the sum scores.  This is useful for regression with an intercept term}
\item{impute}{impute="median" replaces missing values with the item medians, impute = "mean" replaces values with the mean response. impute="none" the subject's scores are based upon the average of the keyed, but non missing scores. impute = "none" is probably more appropriate for a large number of missing cases (e.g., SAPA data).  }
}
\details{Although meant for finding correlation weighted scores using the weights from \code{\link{bestScales}}, it also possible to use alternative weight matrices, such as those returned by the coefficients in \code{\link{lm}}.   
}
\value{
A data frame of scores.}

\author{William Revelle}

\seealso{
 \code{\link{bestScales}} and  \code{\link{setCor}}
}
\examples{
  
#find the weights from a regression model and then apply them to a new set
#derivation of weights from the first 20 cases 
model.lm <- lm(rating ~ complaints + privileges + learning,data=attitude[1:20,])
#or use setCor to find the coefficents
model <- setCor(rating ~ complaints + privileges +learning,data=attitude[1:20,],std=FALSE)
 
 #Apply these to a different set of data (the last 10 cases)
  #note that the regression coefficients need to be a matrix
  scores.lm <- scoreWtd(as.matrix(model.lm$coefficients),attitude[21:30,],sums=TRUE,std=FALSE)
scores <- scoreWtd(model$coefficients,attitude[21:30,],sums=TRUE,std=FALSE)
describe(scores)  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}