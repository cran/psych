\name{bi.bars}
\alias{bi.bars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draw pairs of bargraphs based on two groups}
\description{When showing e.g., age or education distributions for two groups, it is convenient to plot them back to back.  bi.bars will do so.}
\usage{bi.bars(x,grp,horiz,color,...) }

\arguments{
  \item{x}{The data to be drawn}
  \item{grp}{a grouping variable.}
  \item{horiz}{horizontal (default) or vertical bars}
  \item{color}{colors for the two groups -- defaults to blue and red}
  \item{\dots}{Further parameters to pass to the graphing program}
}
\details{A trivial, if useful, function to draw back to back histograms/barplots. One for each group.}
\value{a graphic}

\author{William Revelle}

\examples{
data(bfi)
with(bfi,{bi.bars(age,gender,ylab="Age",main="Age by males and females")
       bi.bars(education,gender,xlab="Education",main="Education by gender",horiz=FALSE)})
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

