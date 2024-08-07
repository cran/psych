\name{bi.bars}
\alias{bi.bars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draw pairs of bargraphs based on two groups}
\description{When showing e.g., age or education distributions for two groups, it is convenient to plot them back to back.  bi.bars will do so.}
\usage{bi.bars(x,var=NULL,grp=NULL,horiz,color,label=NULL,zero=FALSE,xlab,ylab,...) }

\arguments{
  \item{x}{The data frame or matrix from which we specify the data}
  \item{var}{The variable to plot}
  \item{grp}{a grouping variable.}
  \item{horiz}{horizontal (default) or vertical bars}
  \item{color}{colors for the two groups -- defaults to blue and red}
  \item{label}{If specified, labels  for the dependent axis } 
  \item{zero}{If TRUE, subtract the minimum value to make the numbers range from 0 to max -min.  This is useful if  showing heights}
  \item{xlab}{xaxis label if appropriate}
  \item{ylab}{y axis label otherwise}
 \item{\dots}{Further parameters to pass to the graphing program}
}
\details{A trivial, if useful, function to draw back to back histograms/barplots. One for each group.}
\value{a graphic}
\seealso{  \code{\link{describe}}, \code{\link{describeBy}} and  \code{\link{statsBy}} for descriptive statistics and \code{\link{error.bars}} \code{\link{error.bars.by}} and  \code{\link{densityBy}}  \code{\link{violinBy}} for graphic displays }
\author{William Revelle}

\examples{
#data(bfi)
bi.bars(bfi,"age","gender" ,ylab="Age",main="Age by males and females")
 bi.bars(bfi,"education","gender",xlab="Education",
     main="Education by gender",horiz=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

