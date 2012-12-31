\name{withinBetween}
\alias{withinBetween}
\docType{data}
\title{An example of the distinction between within group and between group correlations}
\description{A demonstration that a correlation may be decomposed to a within group correlation  and a between group correlations and these two correlations are independent.  Between group correlations are sometimes called ecological correlations, the decomposition into  within and between group correlations is a basic concept in multilevel modeling.  This data set shows the composite correlations between 9 variables, representing 16 cases with four groups.    
}
\usage{data(withinBetween)}
\format{
  A data frame with 16 observations on the following 10 variables.
  \describe{
 \item{\code{Group}}{An example grouping factor.}
 \item{\code{V1}}{A column of 16 observations }
 \item{\code{V2}}{A column of 16 observations }
 \item{\code{V3}}{A column of 16 observations }
 \item{\code{V4}}{A column of 16 observations }
 \item{\code{V5}}{A column of 16 observations }
 \item{\code{V6}}{A column of 16 observations }
 \item{\code{V7}}{A column of 16 observations }
 \item{\code{V8}}{A column of 16 observations }
 \item{\code{V9}}{A column of 16 observations }
                    
     }
}
\details{Correlations between individuals who belong to different natural groups (based upon e.g., ethnicity, age, gender, college major,or country) reflect an unknown mixture of the pooled correlation within each group as well as the correlation of the means of these groups. These two correlations are independent and do not allow inferences from one level (the group) to the other level (the individual).  This data set shows this independence.  The within group correlations between 9 variables are set to be 1, 0, and -1 while those between groups are also set to be 1, 0, -1.  These two sets of correlations are crossed such that V1, V4, and V7 have within group correlations of 1, as do V2, V5 and V8, and V3, V6 and V9.  V1 has a within group correlation of 0 with V2, V5, and V8, and a -1 within group correlation with V3, V6 and V9.  V1, V2, and V3 share a between group correlation of 1, as do V4, V5 and V6, and V7, V8 and V9.  The first group has a 0 between group correlation with the second and a -1 with the third group.  

\code{\link{statsBy}} can decompose the observed correlation in the between and within correlations.  \code{\link{sim.multilevel}} can produce similar data.

}
\source{The data were created for this example}
\references{
P. D. Bliese. Multilevel modeling in R (2.3) a brief introduction to R, the multilevel package and the nlme package, 2009.

Pedhazur, E.J. (1997) Multiple regression in behavioral research: explanation and prediction.  Harcourt Brace.

Revelle, W. An introduction to psychometric theory with applications in R (in prep) Springer. Draft chapters available at  \url{http://personality-project.org/r/book/}

}

\seealso{ \code{\link{statsBy}},   \code{\link{describeBy}}, and  \code{\link{sim.multilevel}} 
 }
\examples{
data(withinBetween)
pairs.panels(withinBetween,bg=c("red","blue","white","black")[withinBetween[,1]],pch=21,ellipses=FALSE)
stats <- statsBy(withinBetween,'Group')
print(stats,short=FALSE)
}
\keyword{datasets}
