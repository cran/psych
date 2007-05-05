\name{kurtosi}
\alias{kurtosi}

\title{ Kurtosis of a vector, matrix, or data frame}
\description{Find the kurtosis of a vector, matrix, or dataframe.
}
\usage{
kurtosi(x, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector, matrix, or data frame }
  \item{na.rm}{ na.rm =TRUE removes missing data from the column  }
}
\details{Kurtosis in the E1071 program finds the kurtosis for a single vector.  This does it for matrices and dataframes. Used in the describe function
}
\value{
  \item{kurtosi}{a vector of the kurtosis for each column of the matrix}
 
}

\author{William Revelle }

\seealso{  \code{\link{skew}}, \code{\link{describe}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
 round(kurtosi(attitude),2)

## The function is currently defined as
function (x, na.rm = TRUE) 
{
    if (length(dim(x)) == 0) {
        if (na.rm) {
            x <- x[!is.na(x)]
        			}
        
        mx <- mean(x)
         sdx <- sd(x,na.rm=na.rm)
        kurt <- sum((x - mx)^4)/(length(x) * sd(x)^4)  -3
        } else {
    
    kurt <- rep(NA,dim(x)[2])
    mx <- colMeans(x,na.rm=na.rm)
     sdx <- sd(x,na.rm=na.rm)
    for (i in 1:dim(x)[2]) {
    kurt[i] <- sum((x[,i] - mx[i])^4,  na.rm = na.rm)/((length(x[,i]) - sum(is.na(x[,i]))) * sdx[i]^4)  -3
            }
    }
    return(kurt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }

