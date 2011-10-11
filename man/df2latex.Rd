\name{df2latex}
\alias{df2latex}

\title{Convert a data frame or matrix to a LaTeX table}

\description{A handy helper function to convert data frames or matrices to LaTeX tables. Although Sweave is the preferred means of converting R output to LaTeX, it is sometimes useful to go directly from a data.frame or matrix to a LaTeX table.  
}
\usage{
df2latex(x,digits=2,rowlabels=TRUE,caption="df2latex")
}

\arguments{
  \item{x}{A data frame or matrix to convert to LaTeX}
  \item{digits}{Round the output to digits of accuracy.  NULL for formatting character data}
  \item{rowlabels}{If TRUE, use the row names from the matrix or data.frame}
  \item{caption}{The table caption}
  
}
\value{A LaTeX table}

\author{William Revelle
}

\seealso{
The many LaTeX conversion routines in Hmisc.
}
\examples{
df2latex(Thurstone,rowlabels=FALSE)
df2latex(describe(sat.act)[2:10])
 }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }

