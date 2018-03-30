\name{cor.ci}
\alias{cor.ci}
\alias{corCi}
\title{Bootstrapped confidence intervals for raw and composite correlations}

\description{Although normal theory provides confidence intervals for correlations, this is particularly problematic with Synthetic Aperture Personality Assessment (SAPA) data where the individual items are Massively Missing at Random.  Bootstrapped confidence intervals are found for Pearson, Spearman, Kendall, tetrachoric, or polychoric correlations and for scales made from those correlations.
}
\usage{
cor.ci(x, keys = NULL, n.iter = 100,  p = 0.05,overlap = FALSE, 
 poly = FALSE, method = "pearson", plot=TRUE,minlength=5,...)
 corCi(x, keys = NULL, n.iter = 100,  p = 0.05,overlap = FALSE, 
 poly = FALSE, method = "pearson", plot=TRUE,minlength=5,...)

}

\arguments{
  \item{x}{The raw data}
  \item{keys}{If NULL, then the confidence intervals of the raw correlations are found.  Otherwise, composite scales are formed from the keys applied to the correlation matrix (in a logic similar to \code{\link{cluster.cor}} but without the bells and whistles) and the confidence of those composite scales intercorrelations. }
  \item{n.iter}{The number of iterations to bootstrap over. This will be very slow if using tetrachoric/or polychoric correlations. }
  \item{p}{The upper and lower confidence region will include 1-p of the distribution.}
  \item{overlap}{If true, the correlation between overlapping scales is corrected for item overlap.}
  \item{poly}{if FALSE, then find the correlations using the method specified (defaults to Pearson).  If TRUE, the polychoric correlations will be found (slowly).  Because the polychoric function uses multicores (if available), and cor.ci does as well, the number of cores used is options("mc.cores")^2. }
  \item{method}{"pearson","spearman", "kendall"}
  \item{plot}{Show the correlation plot with correlations scaled by the probability values.  To show the matrix in terms of the confidence intervals, use \code{\link{cor.plot.upperLowerCi}}.}
 \item{minlength}{What is the minlength to use in abbreviations of the cis? Defaults to 5}
  \item{...}{Other parameters for axis (e.g., cex.axis to change the font size, srt to rotate the numbers in the plot)}
  }
\details{
The original data are and correlations are found.  If keys are specified (the normal case), then composite scales based upon the correlations are found and reported.  This is the same procedure as done using \code{\link{cluster.cor}} or \code{\link{scoreItems}}.

Then, n.iter times, the data are recreated by sampling subjects (rows) with replacement and the correlations (and composite scales) are found again (and again and again).  Mean and standard deviations of these values are calculated based upon the Fisher Z transform of the correlations.  Summary statistics include the original correlations and their confidence intervals.  For those who want the complete set of replications, those are available as an object in the resulting output.

Although particularly useful for SAPA (\url{https://sapa-project.org}) type data, this will work for any normal data set as well. 

Although the correlations are shown automatically as a \code{\link{cor.plot}}, it is possible to show the upper and lower confidence intervals by using \code{\link{cor.plot.upperLowerCi}}. This will also return, invisibly, a matrix for printing with the lower and upper bounds of the correlations shown below and above the diagonal.
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



\seealso{\code{\link{make.keys}}, \code{\link{cluster.cor}}, and \code{\link{scoreItems}} for forming synthetic correlation matrices from composites of item correlations.  See \code{\link{scoreOverlap}} for correcting for item overlap in scales. See also \code{\link{corr.test}} for standard significance testing of correlation matrices.  See also \code{\link{lowerCor}} for finding and printing correlation matrices, as well as \code{\link{lowerMat}} for displaying them. Also see \code{\link{cor.plot.upperLowerCi}} for displaying the confidence intervals graphically.
}
\examples{
cor.ci(bfi[1:200,1:10])  # just the first 10 variables
#The keys have overlapping scales
  keys.list <- list(agree=c("-A1","A2","A3","A4","A5"), conscientious= c("C1", 
  "C2","C3","-C4","-C5"),extraversion=c("-E1","-E2","E3","E4","E5"), neuroticism= 
  c("N1", "N2", "N3","N4","N5"), openness = c("O1","-O2","O3","O4","-O5"), 
  alpha=c("-A1","A2","A3","A4","A5","C1","C2","C3","-C4","-C5","N1","N2","N3","N4","N5"),
beta = c("-E1","-E2","E3","E4","E5","O1","-O2","O3","O4","-O5") )
  keys <- make.keys(bfi,keys.list)
  
#do not correct for item overlap
rci <-  cor.ci(bfi[1:200,],keys,n.iter=10,main="correlation with overlapping scales") 
#also shows the graphic -note the overlap
#correct for overlap
rci <-  cor.ci(bfi[1:200,],keys,overlap=TRUE, n.iter=10,main="Correct for overlap") 
#show the confidence intervals
ci <- cor.plot.upperLowerCi(rci)  #to show the upper and lower confidence intervals
ci   #print the confidence intervals in matrix form
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{ models }