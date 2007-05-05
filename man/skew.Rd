\name{skew}
\alias{skew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate skew for a vector, matrix, or data.frame}
\description{Find the skew for each variable in a data.frame or matrix.  Unlike skew in e1071, this reports a different skew for each variable.
}
\usage{
skew(x, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data.frame or matrix }
  \item{na.rm}{ how to treat missing data }
}
\details{given a matrix or data.frame x, find the skew for each column.
}
\value{if input is a matrix or data.frame, skew is a vector of skews}

\author{ William Revelle}

\seealso{\code{\link{describe}} }
\examples{
#round(skew(attitude),2)
## The function is currently defined as
function (x, na.rm = TRUE) {
    if(length(dim(x))==0) {
       if (na.rm)   { x <- x[!is.na(x)] }        #remove missing values  
    sum((x - mean(x))^3)/(length(x) * sd(x)^3)  #calculate skew  for a vector
 } else { apply(x,2,function(x) sum((x - mean(x,na.rm=na.rm))^3,na.rm=na.rm)/( (length(x)-sum(is.na(x)) )* sd(x,na.rm=na.rm)^3)) }  }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
