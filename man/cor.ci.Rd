\name{cor.ci}
\alias{cor.ci}
\title{Bootstrapped confidence intervals for raw and composite correlations}

\description{Although normal theory provides confidence intervals for correlations, this is particularly problematic with Synthetic Aperture Personality Assessment (SAPA) data where the individual items are Massively Missing at Random.  Bootstrapped confidence intervals are found for Pearson, Spearman, Kendall, tetrachoric, or polychoric correlations and for scales made from those correlations.
}
\usage{
cor.ci(x, keys = NULL, n.iter = 100,  p = 0.05, poly = FALSE, method = "pearson")
}

\arguments{
  \item{x}{The raw data}
  \item{keys}{If NULL, then the confidence intervals of the raw correlations are found.  Otherwise, composite scales are formed from the keys applied to the correlation matrix (in a logic similar to \code{\link{cluster.cor}} but without the bells and whistles) and the confidence of those composite scales intercorrelations. }
  \item{n.iter}{The number of iterations to bootstrap over. This will be very slow if using tetrachoric/or polychoric correlations. }
  \item{p}{The upper and lower confidence region will include 1-p of the distribution.}
  \item{poly}{if FALSE, then find the correlations using the method specified (defaults to Pearson).  If TRUE, the polychoric correlations will be found (slowly)}
  \item{method}{"pearson","spearman", "kendall"}
  }
\details{
The original data are and correlations are found.  If keys are specified (the normal case), then composite scales based upon the correlations are found and reported.  This is the same procedure as done using \code{\link{cluster.cor}} or \code{\link{score.items}}.

Then, n.iter times, the data are recreated by sampling subjects (rows) with replacement and the correlations (and composite scales) are found again (and again and again).  Mean and standard deviations of these values are calculated based upon the Fisher Z transform of the correlations.  Summary statistics include the original correlations and their confidence intervals.  For those who want the complete set of replications, those are available as an object in the resulting output.

Although particularly useful for SAPA (\url{http://sapa-project.org}) type data, this will work for any normal data set as well. 
}
\value{

 \item{rho }{The original (composite) correlation matrix. }
\item{means }{Mean (Fisher transformed) correlation}
\item{sds}{Standard deviation of Fisher transformed correlations}
\item{ci}{Mean +/- alpha/2 of the z scores as well as the alpha/2 and 1-alpha/2 quantiles. These are labeled as lower.emp(ircal), lower.norm(al), upper.norm and upper.emp.}
\item{replicates}{The observed replication values so one can do one's own estimates}
}
\references{For SAPA type data, see Revelle, W., Wilt, J.,  and Rosenthal, A. (2010)  Personality and Cognition: The Personality-Cognition Link. In Gruszka, A.  and Matthews, G. and Szymura, B. (Eds.) Handbook of Individual Differences in Cognition: Attention, Memory and Executive Control, Springer. 
}
\author{William Revelle}



\seealso{\code{\link{make.keys}}, \code{\link{cluster.cor}}, and \code{\link{score.items}} for forming synthetic correlation matrices from composites of item correlations.  See also \code{\link{corr.test}} for standard significance testing of correlation matrices.  See also \code{\link{lowerCor}} for finding and printing correlation matrices, as well as \code{\link{lowerMat}} for displaying them. 
}
\examples{
 keys.list <-
  list(agree=c("-A1","A2","A3","A4","A5"),conscientious=c("C1","C2","C2","-C4","-C5"),
extraversion=c("-E1","-E2","E3","E4","E5")) 
  keys <- make.keys(bfi[1:15],keys.list)
  cor.ci(bfi[1:15],keys,n.iter=10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{ models }