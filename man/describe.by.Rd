\name{describeBy}
\alias{describeBy}
\alias{describe.by}

\title{ Basic summary statistics by group}
\description{Report basic summary statistics by a grouping variable.  Useful if the grouping variable is some experimental variable and data are to be aggregated for plotting.  Partly a wrapper for by and \code{\link{describe}}
}
\usage{
describeBy(x, group=NULL,mat=FALSE,type=3,digits=15,...)
describe.by(x, group=NULL,mat=FALSE,type=3,...)  # deprecated

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data.frame or matrix. See note for statsBy.  }
  
  \item{group}{a grouping variable or a list of grouping variables}
  \item{mat}{provide a matrix output rather than a list}
  \item{type}{Which type of skew and kurtosis should be found}
  \item{digits}{When giving matrix output, how many digits should be reported?}
  
  \item{...}{parameters to be passed to describe}
}
\details{To get descriptive statistics for several different grouping variables, make sure that group is a list.  In the case of matrix output with multiple grouping variables, the grouping variable values are added to the output.

The type parameter specifies which version of skew and kurtosis should be found.  See \code{\link{describe}} for more details. 

An alternative function (\code{\link{statsBy}}) returns a list of means, n, and standard deviations for each group.  This is particularly useful if finding weighted correlations of group means using \code{\link{cor.wt}}. More importantly, it does a proper within and between group decomposition of the correlation.

\code{\link{cohen.d}}  will work for two groups. It converts the data into mean differences and pools the within group standard deviations.  Returns cohen.d statistic as well as the multivariate generalization (Mahalanobis D). 
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



\seealso{  \code{\link{describe}},  \code{\link{statsBy}}, \code{\link{densityBy}} and \code{\link{violinBy}},  \code{\link{cohen.d}},  \code{\link{cohen.d.by}}, and  \code{\link{cohen.d.ci}}
 as well as \code{\link{error.bars}} and \code{\link{error.bars.by}} for other graphical displays. }
\examples{

data(sat.act)
describeBy(sat.act,sat.act$gender) #just one grouping variable	
#describeBy(sat.act,list(sat.act$gender,sat.act$education))  #two grouping variables
des.mat <- describeBy(sat.act$age,sat.act$education,mat=TRUE) #matrix (data.frame) output 
des.mat <- describeBy(sat.act$age,list(sat.act$education,sat.act$gender),
             mat=TRUE,digits=2)  #matrix output


}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ univar }% __ONLY ONE__ keyword per line
