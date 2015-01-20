\name{mediate}
\alias{mediate}
\alias{mediate.diagram}
\title{Estimate and display direct and indirect effects of mediators and  moderator in path models}
\description{
Find the direct and indirect effects of a predictor in path models of mediation and moderation. Bootstrap confidence intervals for the indirect effects.  
}
\usage{
mediate(y, x, m, data, mod = NULL, n.obs = NULL, use = "pairwise", n.iter = 5000, 
       alpha = 0.05, std = FALSE)
mediate.diagram(medi,digits=2,...)
}

\arguments{
  \item{y}{The dependent variable (or a formula suitable for a linear model)}
  \item{x}{The principal Independent Variable}
 \item{m}{One (or more) mediating variables}
  \item{data}{A data frame holding the data or a correlation  matrix. }
  \item{mod}{A moderating variable, if desired}
  \item{n.obs}{If the data are from a correlation or covariance matrix, 
     how many observations were used.}
  \item{use}{use="pairwise" is the default when finding correlations or covariances}
  \item{n.iter}{Number of bootstrap resamplings to conduct}
  \item{alpha}{Set the width of the confidence interval to be 1 - alpha}
  \item{std}{standardize the covariances to find the standardized betas}
  \item{digits}{The number of digits to report in the mediate.diagram.} 
  \item{medi}{The output from mediate may be imported into mediate.diagram}
  \item{...}{Additional graphical parameters to pass to mediate.diagram}
}
\details{
When doing linear modeling, it is frequently convenient to estimate the direct effect of a predictor controlling for the indirect effect of a mediator.  See Preacher and Hayes (2004) for a very thorough discussion of mediation
}
\value{
\item{total}{The total direct effect of x on y  (c)}
\item{direct}{The beta effects of x (c') and m  (b) on y   }
\item{indirect}{The indirect effect of x through m on y (c-ab)}
\item{mean.boot}{mean bootstrapped value of indirect effect}
\item{sd.boot}{Standard deviation of bootstrapped values}
\item{ci.quant}{The upper and lower confidence intervals based upon the quantiles of the bootstrapped distribution.}
\item{boot}{The bootstrapped values themselves.}
\item{a}{The effect of x on m}
\item{b}{The effect of m on y}
\item{b.int}{The interaction of x and mod (if specified)}
}
\references{
Preacher, Kristopher J and Hayes, Andrew F (2004) SPSS and SAS procedures for estimating indirect effects in simple mediation models. Behavior Research Methods, Instruments, \& Computers 36, (4) 717-731.

Data from Preacher and Hayes (2004).

}
\author{William Revelle
}
\note{Beta version.  Currently just does simple (and parallel) mediation and simple moderation.  The graphics are also fairly limited. 
}


\seealso{\code{\link{setCor}} and \code{\link{setCor.diagram}}

}
\examples{
#data from Preacher and Hayes (2004)
sobel <- structure(list(SATIS = c(-0.59, 1.3, 0.02, 0.01, 0.79, -0.35, 
-0.03, 1.75, -0.8, -1.2, -1.27, 0.7, -1.59, 0.68, -0.39, 1.33, 
-1.59, 1.34, 0.1, 0.05, 0.66, 0.56, 0.85, 0.88, 0.14, -0.72, 
0.84, -1.13, -0.13, 0.2), THERAPY = structure(c(0, 1, 1, 0, 1, 
1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 
1, 1, 1, 0), value.labels = structure(c(1, 0), .Names = c("cognitive", 
"standard"))), ATTRIB = c(-1.17, 0.04, 0.58, -0.23, 0.62, -0.26, 
-0.28, 0.52, 0.34, -0.09, -1.09, 1.05, -1.84, -0.95, 0.15, 0.07, 
-0.1, 2.35, 0.75, 0.49, 0.67, 1.21, 0.31, 1.97, -0.94, 0.11, 
-0.54, -0.23, 0.05, -1.07)), .Names = c("SATIS", "THERAPY", "ATTRIB"
), row.names = c(NA, -30L), class = "data.frame", variable.labels = structure(c("Satisfaction", 
"Therapy", "Attributional Positivity"), .Names = c("SATIS", "THERAPY", 
"ATTRIB")))

mediate(1,2,3,sobel)  #The example in Preacher and Hayes

}

\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
