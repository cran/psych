\name{structure.list}
\alias{structure.list}
\alias{phi.list}

\title{Create factor model matrices from an input list}
\description{When creating a structural diagram or a structural model, it is convenient to not have to specify all of the zero loadings in a structural matrix.  structure.list converts list input into a design matrix.  phi.list does the same for a correlation matrix. Factors with NULL values are filled with 0s.
}
\usage{
structure.list(nvars, f.list,f=NULL, f.labels = NULL, item.labels = NULL)
phi.list(nf,f.list, f.labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvars}{Number of variables in the design matrix }
  \item{f.list}{A list of items included in each factor (for structure.list, or the factors that correlate with the specified factor for phi.list}
  \item{f}{prefix for parameters -- needed in case of creating an X set and a Y set}
  \item{f.labels}{Names for the factors }
  \item{item.labels}{Item labels }
  
  \item{nf}{Number of factors in the phi matrix}
}
\details{This is almost self explanatory.  See the examples.
}
\value{
  \item{factor.matrix}{a matrix of factor loadings to model}
}

\seealso{  \code{\link{structure.graph}}}
\examples{
fx <- structure.list(9,list(F1=c(1,2,3),F2=c(4,5,6),F3=c(7,8,9)))
fy <- structure.list(3,list(Y=c(1,2,3)),"Y")
phi <- phi.list(4,list(F1=c(4),F2=c(1,4),F3=c(2),F4=c(1,2,3)))
fx
phi
fy


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{models }