\name{alpha}
\alias{alpha}
\alias{alpha.scale}

\title{Find two estimates of reliability: Cronbach's alpha and Guttman's Lambda 6. }
\description{
 Internal consistency measures of reliability range from \eqn{\omega_h}{omega_hierchical} to \eqn{\alpha}{alpha} to \eqn{\omega_t}{omega_total}.  This function reports two estimates: Cronbach's coefficient \eqn{\alpha}{alpha} and Guttman's \eqn{\lambda_6}{lambda_6}.  Also reported are item - whole correlations, \eqn{\alpha}{alpha} if an item is omitted, and item means and standard deviations.
}
\usage{
alpha(x, keys=NULL,cumulative=FALSE, title=NULL, na.rm = TRUE)
}

\arguments{
  \item{x}{A data.frame or matrix of data, or a covariance or correlation matrix }
  \item{keys}{If some items are to be reversed keyed, then the direction of all items must be specified in a keys vector}
  \item{title}{Any text string to identify this run}
  \item{cumulative}{should means reflect the sum of items or the mean of the items.  The default value is means.}
  \item{na.rm}{The default is to remove missing values and find pairwise correlations}
}
\details{Alpha is one of several estimates of the internal consistency reliability of a test.

Surprisingly, 105 years after Spearman (1904) introduced the concept of reliability to psychologists, there are still multiple approaches for measuring it. Although very popular, Cronbach's \eqn{\alpha}  (1951) underestimates the reliability of a test and over estimates the first factor saturation.

\eqn{\alpha}{alpha} (Cronbach, 1951) is the same as Guttman's  \eqn{\lambda}{lambda}3 (Guttman, 1945) and may be found by
\deqn{
\lambda_3 =  \frac{n}{n-1}\Bigl(1 - \frac{tr(\vec{V})_x}{V_x}\Bigr) = \frac{n}{n-1} \frac{V_x - tr(\vec{V}_x)}{V_x} = \alpha
}{Lambda 3 = (n)/(n-1)(1-tr(Vx)/(Vx)  = (n)/(n-1)(Vx-tr(Vx)/Vx  = alpha}

Perhaps because it is so easy to calculate and is available in most commercial programs, alpha is without doubt the most frequently reported measure of internal consistency reliability. Alpha is the mean of all possible spit half reliabilities (corrected for test length).  For a unifactorial test, it is a reasonable estimate of the first factor saturation, although if the test has any microstructure (i.e., if it is ``lumpy") coefficients \eqn{\beta}{beta} (Revelle, 1979; see \code{\link{ICLUST}}) and \eqn{\omega_h}{omega_hierchical} (see \code{\link{omega}})  are more appropriate estimates of the general factor saturation.  \eqn{\omega_t}{omega_total} (see \code{\link{omega}}) is a better estimate of the reliability of the total test.  

Guttman's Lambda 6 (G6) considers the amount of variance in each item that can be accounted for the linear regression of all of the other items (the squared multiple correlation or smc), or more precisely, the variance of the errors, \eqn{e_j^2},  and is
\deqn{
\lambda_6 = 1 - \frac{\sum e_j^2}{V_x} = 1 - \frac{\sum(1-r_{smc}^2)}{V_x}
.}{lamada 6 = 1 - sum(e^2)/Vx = 1-sum(1-r^2(smc))/Vx.}

The squared multiple correlation is a lower bound for the item communality and as the number of items increases, becomes a better estimate.

G6 is also sensitive to lumpyness in the test and should not be taken as a measure of unifactorial structure.  For lumpy tests, it will be greater than alpha.  For tests with equal item loadings, alpha > G6, but if the loadings are unequal or if there is a general factor, G6 > alpha. 

Alpha and G6 are both positive functions of the number of items in a test as well as the average intercorrelation of the items in the test.  When calculated from the item variances and total test variance, as is done here, raw alpha is sensitive to differences in the item variances. Standardized alpha is based upon the correlations rather than the covariances.  

More complete reliability analyses of a single scale can be done using the \code{\link{omega}} function which finds \eqn{\omega_h}{omega_hierchical} and \eqn{\omega_t}{omega_total} based upon a hierarchical factor analysis.  

Alternative functions \code{\link{score.items}} and  \code{\link{cluster.cor}} will also score multiple scales and report more useful statistics. ``Standardized" alpha is calculated from the inter-item correlations and will differ from raw alpha. 

}
\value{
  \item{total }{a list containing}
  \item{raw_alpha}{alpha based upon the covariances}
  \item{std.alpha}{The standarized alpha based upon the correlations}
  \item{G6(smc)}{Guttman's Lambda 6 reliability}
  \item{average_r}{The average interitem correlation}
  \item{mean}{For data matrices, the mean of the scale formed by summing the items}
  \item{sd}{For data matrices, the standard deviation of the total score}
  \item{alpha.drop }{A data frame with all of the above for the case of each item being removed one by one.}
  \item{item.stats}{A data frame including}
  \item{r}{The correlation of each item with the total score (not corrected for item overlap)}
  \item{r.cor}{Item whole correlation corrected for item overlap and scale reliability}
  \item{mean}{for data matrices, the mean of each item}
  \item{sd}{For data matrices, the standard deviation of each item}
  
}
\references{
Cronbach, L.J. (1951) Coefficient alpha and the internal strucuture of tests.  Psychometrika, 16, 297-334.

Guttman, L. (1945). A basis for analyzing test-retest reliability. Psychometrika, 10 (4), 255-282. 

Revelle, W. Hierarchical Cluster Analysis and the Internal Structure of Tests. Multivariate Behavioral Research, 1979, 14, 57-74.

Revelle, W. and Zinbarg, R. E. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma. Psychometrika, 2009. 


	
}
\author{William Revelle }


\seealso{  \code{\link{omega}}, \code{\link{ICLUST}}, \code{\link{guttman}}, \code{\link{score.items}}, \code{\link{cluster.cor}}}
\examples{
r4 <- sim.congeneric()
alpha(r4)
r9 <- sim.hierarchical()
alpha(r9)
#an example of two independent factors that produce reasonable alphas
#this is a case where alpha is a poor indicator of unidimensionality
two.f <- sim.item(8)
alpha(two.f,keys=c(rep(1,4),rep(-1,4))) 


}

\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line
