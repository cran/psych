\name{skew}
\alias{skew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate skew for a vector, matrix, or data.frame}
\description{Find the skew for each variable in a data.frame or matrix.  Unlike skew in e1071, this calculates a different skew for each variable.
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

\seealso{\code{\link{describe}}, \code{\link{describe.by}}, \code{\link{kurtosi}}   }
\examples{
round(skew(attitude),2)
## The function is currently defined as
function (x, na.rm = TRUE) 
{
    if (length(dim(x)) == 0) {
        if (na.rm) {
            x <- x[!is.na(x)]
        			}
        sdx <- sd(x,na.rm=na.rm)
        mx <- mean(x)
        skewer <- sum((x - mx)^3)/(length(x) * sd(x)^3) 
        } else {
    
    skewer <- rep(NA,dim(x)[2])
    mx <- colMeans(x,na.rm=na.rm)
    sdx <- sd(x,na.rm=na.rm)
    for (i in 1:dim(x)[2]) {
    skewer[i] <- sum((x[,i] - mx[i])^3,  na.rm = na.rm)/((length(x[,i]) - sum(is.na(x[,i]))) * sdx[i]^3)
            }
    }
    return(skewer)
}
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
