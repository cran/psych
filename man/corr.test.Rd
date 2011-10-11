\name{corr.test}
\alias{corr.test}
\alias{corr.p}
\title{Find the correlations, sample sizes, and probability values between elements of a matrix or data.frame.   }
\description{Although the cor function finds the correlations for a matrix, it does not report probability values. corr.test uses cor to find the correlations for either complete or pairwise data and reports the sample sizes and probability values as well. For symmetric matrices, raw probabilites are reported below the diagonal and correlations adjusted for multiple comparisons above the diagonal.
}
\usage{
corr.test(x, y = NULL, use = "pairwise",method="pearson",adjust="holm")
corr.p(r,n,adjust="holm")
}

\arguments{
  \item{x}{A matrix or dataframe }
  \item{y}{A second matrix or dataframe with the same number of rows as x }
  \item{use}{use="pairwise" is the default value and will do pairwise deletion of cases. use="complete" will select just complete cases. }
  \item{method}{method="pearson" is the default value.  The alternatives to be passed to cor are "spearman" and "kendall"}
  \item{adjust}{What adjustment for multiple tests should be used? ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). See \code{\link{p.adjust}} for details about why to use "holm" rather than "bonferroni"). }
  \item{r}{A correlation matrix}
  \item{n}{Number of observations of using corr.p}
}
\details{corr.test uses the \code{\link{cor}} function to find the correlations, and then applies a t-test to the individual correlations using the formula
\deqn{t = \frac{r * \sqrt(n-2)}{\sqrt(1-r^2)}
}{t = r* sqrt(n-1)/sqrt(1-r^2) }

The probability values may be adjusted using the Holm (or other) correction.  If the matrix is symmetric (no y data), then the original p values are reported below the diagonal and the adjusted above the diagonal.  Otherwise, all probabilities are adjusted (unless adjust="none").
}
\value{
  \item{r}{The matrix of correlations}
  \item{n}{Number of cases per correlation}
  \item{t}{value of t-test for each correlation}
  \item{p}{two tailed probability of t for each correlation.  For symmetric matrices, p values adjusted for multiple tests are reported above the diagonal. }
}

\seealso{  \code{\link{cor.test}} for tests of a single correlation,  Hmisc::rcorr for an equivalant function, \code{\link{r.test}} to test the difference between correlations, and \code{\link{cortest.mat}} to test for equality of two correlation matrices. 

In particular, see \code{\link{p.adjust}} for a discussion of p values associated with multiple tests.}
\examples{
data(sat.act)
corr.test(sat.act)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{ models }% __ONLY ONE__ keyword per line
