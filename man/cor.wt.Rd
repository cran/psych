\name{cor.wt}
\alias{cor.wt}
\title{The sample size weighted correlation may be used in correlating aggregated data}
\description{If using aggregated data, the correlation of the means does not reflect the sample size used for each mean. cov.wt in RCore does this and returns a covariance matrix or the correlation matrix.  The cor.wt function weights by sample size or by standard errors and by default return correlations. 
}
\usage{
cor.wt(data,vars=NULL, w=NULL,sds=NULL, cor=TRUE)
}

\arguments{
  \item{data}{A matrix or data frame}
  \item{vars}{Variables to analyze}
 \item{w}{A set of weights (e.g., the sample sizes)}
  \item{sds}{Standard deviations of the samples (used if weighting by standard errors)}
  \item{cor}{Report correlations (the default) or covariances}
}
\details{A weighted correlation is just \eqn{ r_{ij} = \frac{\sum(wt_{k} (x_{ik} - x_{jk})}{\sqrt{wt_{ik} \sum(x_{ik}^2) wt_jk \sum(x_{jk}^2)}}    }{\sum (wt_k * (x_ik - x_jk)) /sqrt[wt_k \sum(x^2_ik) wt_k \sum(x^2_jk)]}  where \eqn{x_{ik}}{x_ik} is a deviation from the weighted mean.  

The weighted correlation is appropriate for correlating aggregated data, where individual data points might reflect the means of a number of observations.  In this case, each point is weighted by its sample size (or alternatively, by the standard error).  If the weights are all equal, the correlation is just a normal Pearson correlation. 

Used when finding correlations of group means found using \code{\link{statsBy}}.
} 
\value{
   \item{cor }{The weighted correlation}
  \item{xwt}{The data as weighted deviations from the weighted mean }
  \item{wt}{The weights used (calculated from the sample sizes).}
  \item{mean}{The weighted means}
  \item{xc}{Unweighted, centered deviation scores from the weighted mean}
  \item{xs}{Deviation scores weighted by the standard error of each sample mean}
  }

\author{William Revelle
}
\note{A generalization of \code{\link{cov.wt}} in core R}


\seealso{
 See Also as \code{\link{cov.wt}}, \code{\link{statsBy}}
}
\examples{
means.by.age <- statsBy(sat.act,"age")
wt.cors <- cor.wt(means.by.age)
lowerMat(wt.cors$r)  #show the weighted correlations
unwt <- lowerCor(means.by.age$mean)
mixed <- lowerUpper(unwt,wt.cors$r)  #combine both results
cor.plot(mixed,TRUE,main="weighted versus unweighted correlations")
diff <- lowerUpper(unwt,wt.cors$r,TRUE)
cor.plot(diff,TRUE,main="differences of weighted versus unweighted correlations")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
