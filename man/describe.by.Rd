\name{describe.by}
\alias{describe.by}
\title{ Basic summary statistics by group}
\description{Report basic summary statistics by a grouping variable.  Useful if the grouping variable is some experimental variable and data are to be aggregated for plotting.  Partly a wrapper for by and \code{\link{describe}}
}
\usage{
describe.by(x, group=NULL,mat=FALSE,type=3,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data.frame or matrix  }
  \item{group}{a grouping variable or a list of grouping variables}
  \item{mat}{provide a matrix output rather than a list}
  \item{type}{Which type of skew and kurtosis should be found}
  \item{...}{parameters to be passed to describe}
}
\details{To get descriptive statistics for several different grouping variables, make sure that group is a list.  In the case of matrix output with multiple grouping variables, the grouping variable values are added to the output.

The type parameter specifies which version of skew and kurtosis should be found.  See \code{\link{describe}} for more details. 
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
des.mat <- describe.by(sat.act$age,sat.act$education,mat=TRUE) #matrix  (data.frame) output 
des.mat <- describe.by(sat.act$age,list(sat.act$education,sat.act$gender),mat=TRUE)  #matrix output


}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ univar }% __ONLY ONE__ keyword per line
