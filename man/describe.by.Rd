\name{describe.by}
\alias{describe.by}
\title{ Basic summary statistics by group}
\description{Report basic summary statistics by a grouping variable.  Useful if the grouping variable is some experimental variable and data are to be aggregated for plotting.  Just a wrapper for by and \code{\link{describe}}.
}
\usage{
describe.by(x, group,mat=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data.frame or matrix  }
  \item{group}{a grouping variable or a list of grouping variables}
  \item{mat}{provide a matrix output rather than a list}
  \item{...}{parameters to be passed to describe}
}
\details{To get descriptive statistics for several different grouping variables, make sure that group is a list!
}
\value{
 A data.frame of the relevant statistics broken down by group: \cr
   item name \cr
   item number \cr
   number of valid cases\cr
   mean\cr
   standard deviation\cr
   median\cr
   mad: median absolute deviation (from the median) \cr
   minimum\cr
   maximum\cr
   skew\cr
   standard error\cr
}

\author{ William Revelle}

\seealso{  \code{\link{describe}}}
\examples{
data(sat.act)
describe.by(sat.act,sat.act$gender) #just one grouping variable	
#describe.by(sat.act,list(sat.act$gender,sat.act$education))  #two grouping variables
#des.mat <- describe.by(sat.act$age,sat.act$education,mat=TRUE) #matrix output 
#des.mat <- describe.by(sat.act$age,list(sat.act$education,sat.act$gender),mat=TRUE) 


}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ univar }% __ONLY ONE__ keyword per line
