\name{df2latex}
\alias{df2latex}
\alias{cor2latex}
\alias{fa2latex}
\alias{omega2latex}
\alias{irt2latex}
\alias{ICC2latex}

\title{Convert a data frame, correlation matrix, or factor analysis output to a LaTeX table}

\description{A set of handy helper functions to convert data frames or matrices to LaTeX tables. Although Sweave is the preferred means of converting R output to LaTeX, it is sometimes useful to go directly from a data.frame or matrix to a LaTeX table.   cor2latex will find the correlations and then create a lower (or upper) triangular matrix for latex output.  fa2latex will create the latex commands for showing the loadings and factor intercorrelations. As the default option, tables are prepared in an approximation of APA format. 
}
\usage{
df2latex(x,digits=2,rowlabels=TRUE,apa=TRUE,short.names=TRUE,font.size ="scriptsize",
big.mark=NULL,drop.na=TRUE, heading="A table from the psych package in R",
caption="df2latex",label="default", char=FALSE,stars=FALSE)
cor2latex(x,use = "pairwise", method="pearson", adjust="holm",stars=FALSE,
       digits=2,rowlabels=TRUE,lower=TRUE,apa=TRUE,short.names=TRUE,
     font.size ="scriptsize",
      heading="A correlation table from the psych package in R.",
      caption="cor2latex",label="default")
fa2latex(f,digits=2,rowlabels=TRUE,apa=TRUE,short.names=FALSE,cumvar=FALSE,
       cut=0,alpha=.05,font.size ="scriptsize",
       heading="A factor analysis table from the psych package in R",
       caption="fa2latex",label="default") 
omega2latex(f,digits=2,rowlabels=TRUE,apa=TRUE,short.names=FALSE,cumvar=FALSE,cut=.2,
        font.size ="scriptsize", 
        heading="An omega analysis table from the psych package in R",
        caption="omega2latex",label="default") 

irt2latex(f,digits=2,rowlabels=TRUE,apa=TRUE,short.names=FALSE,
       font.size ="scriptsize",
     heading="An IRT factor analysis table from R",caption="fa2latex",label="default")
ICC2latex(icc,digits=2,rowlabels=TRUE,apa=TRUE,ci=TRUE,
font.size ="scriptsize",big.mark=NULL, drop.na=TRUE,
    heading="A table from the psych package in R",
   caption="ICC2latex",label="default",char=FALSE) 
}

\arguments{
  \item{x}{A data frame or matrix to convert to LaTeX. If non-square, then correlations will be found prior to printing in cor2latex}
  \item{digits}{Round the output to digits of accuracy.  NULL for formatting character data}
  \item{rowlabels}{If TRUE, use the row names from the matrix or data.frame}
  \item{short.names}{Name the columns 1 .. N to save space}
  \item{apa}{If TRUE formats table in APA style}
  \item{cumvar}{For factor analyses, should we show the cumulative variance accounted for?}
  \item{font.size}{e.g., "scriptsize", "tiny" or anyother acceptable LaTeX font size.}
  \item{heading}{The label appearing at the top of the table}
  \item{caption}{The table caption}
  \item{lower}{in cor2latex, just show the lower triangular matrix}
  \item{f}{The object returned from a factor analysis using \code{\link{fa}} or \code{\link{irt.fa}}. }
  \item{label}{The label for the table}
  \item{big.mark}{Comma separate numbers large numbers (big.mark=",")}
  \item{drop.na}{Do not print NA values}
  \item{method}{When finding correlations, which method should be used (pearson)}
  \item{use}{use="pairwise" is the default when finding correlations in cor2latex}
  \item{adjust}{If showing probabilities, which adjustment should be used (holm)}
  \item{stars}{Should probability asterixs be displayed in cor2latex (FALSE)}
  \item{char}{char=TRUE allows printing tables with character information, but does not allow for putting in commas into numbers}
  \item{cut}{In omega2latex and fa2latex, do not print abs(values) < cut }
  \item{alpha}{If fa has returned confidence intervals, then what values of loadings should be boldfaced?}
  \item{icc}{Either the output of an ICC, or the data to be analyzed.}
  \item{ci}{Should confidence intervals of the ICC be displayed}
}
\value{A LaTeX table.  Note that if showing "stars" for correlations, then one needs to use the siunitx  package in LaTex. }

\author{William Revelle with suggestions from Jason French and David Condon and Davide Morselli}

\seealso{
The many LaTeX conversion routines in Hmisc.
}
\examples{
df2latex(Thurstone,rowlabels=FALSE,apa=FALSE,short.names=FALSE,
        caption="Thurstone Correlation matrix")
df2latex(Thurstone,heading="Thurstone Correlation matrix in APA style")

df2latex(describe(sat.act)[2:10],short.names=FALSE)
cor2latex(Thurstone)
cor2latex(sat.act,short.names=FALSE)
fa2latex(fa(Thurstone,3),heading="Factor analysis from R in quasi APA style")
 }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }

