\name{alpha}
\alias{alpha}
\alias{alpha.scale}
\alias{alpha.ci}

\title{Find two estimates of reliability: Cronbach's alpha and Guttman's Lambda 6. }
\description{
 Internal consistency measures of reliability range from \eqn{\omega_h}{omega_hierchical} to \eqn{\alpha}{alpha} to \eqn{\omega_t}{omega_total}.  This function reports two estimates: Cronbach's coefficient \eqn{\alpha}{alpha} and Guttman's \eqn{\lambda_6}{lambda_6}.  Also reported are item - whole correlations, \eqn{\alpha}{alpha} if an item is omitted, and item means and standard deviations.
}
\usage{
alpha(x, keys=NULL,cumulative=FALSE, title=NULL, max=10,na.rm = TRUE,
   check.keys=FALSE,n.iter=1,delete=TRUE,use="pairwise",warnings=TRUE,
   n.obs=NULL,impute=NULL
   )
alpha.ci(alpha,n.obs,n.var=NULL,p.val=.05,digits=2)
}

\arguments{
  \item{x}{A data.frame or matrix of data, or a covariance or correlation matrix }
  \item{keys}{If some items are to be reversed keyed, then either specify the direction of all items or just a vector of which items to reverse }
  \item{title}{Any text string to identify this run}
  \item{cumulative}{should means reflect the sum of items or the mean of the items.  The default value is means.}
  \item{max}{the number of categories/item to consider if reporting category frequencies.  Defaults to 10, passed to \code{link{response.frequencies}} }
  \item{na.rm}{The default is to remove missing values and find pairwise correlations}
  \item{check.keys}{if TRUE, then find the first principal component and reverse key items with negative loadings.  Give a warning if this happens. }
 \item{n.iter}{Number of iterations if bootstrapped confidence intervals are desired}
  \item{delete}{Delete items with no variance and issue a warning}
  \item{use}{Options to pass to the cor function: "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". The default is "pairwise"}
  \item{warnings}{By default print a warning and a message that items were reversed. Suppress the message if warnings = FALSE}
  \item{alpha}{The value to use for confidence intervals}
  \item{n.obs}{If using correlation matrices as input, by specify the number of observations, we can find confidence intervals}
  \item{impute}{How should we impute missing data? Not at all, medians, or means}
  \item{n.var}{Number of items in the scale (to find r.bar)}
  \item{p.val}{width of confidence interval (pval/2 to 1-p.val/2) }
  \item{digits}{How many digits to use for alpha.ci}
}
\details{Alpha is one of several estimates of the internal consistency reliability of a test.

Surprisingly, more than a century after Spearman (1904) introduced the concept of reliability to psychologists, there are still multiple approaches for measuring it. Although very popular, Cronbach's \eqn{\alpha}  (1951) underestimates the reliability of a test and over estimates the first factor saturation.

\eqn{\alpha}{alpha} (Cronbach, 1951) is the same as Guttman's  \eqn{\lambda}{lambda}3 (Guttman, 1945) and may be found by
\deqn{
\lambda_3 =  \frac{n}{n-1}\Bigl(1 - \frac{tr(\vec{V})_x}{V_x}\Bigr) = \frac{n}{n-1} \frac{V_x - tr(\vec{V}_x)}{V_x} = \alpha
}{Lambda 3 = (n)/(n-1)(1-tr(Vx)/(Vx)  = (n)/(n-1)(Vx-tr(Vx)/Vx  = alpha}

Perhaps because it is so easy to calculate and is available in most commercial programs, alpha is without doubt the most frequently reported measure of internal consistency reliability. Alpha is the mean of all possible spit half reliabilities (corrected for test length).  For a unifactorial test, it is a reasonable estimate of the first factor saturation, although if the test has any microstructure (i.e., if it is ``lumpy") coefficients \eqn{\beta}{beta} (Revelle, 1979; see \code{\link{ICLUST}}) and \eqn{\omega_h}{omega_hierchical} (see \code{\link{omega}})  are more appropriate estimates of the general factor saturation.  \eqn{\omega_t}{omega_total} (see \code{\link{omega}}) is a better estimate of the reliability of the total test.  

Guttman's Lambda 6 (G6) considers the amount of variance in each item that can be accounted for the linear regression of all of the other items (the squared multiple correlation or smc), or more precisely, the variance of the errors, \eqn{e_j^2},  and is
\deqn{
\lambda_6 = 1 - \frac{\sum e_j^2}{V_x} = 1 - \frac{\sum(1-r_{smc}^2)}{V_x}
.}{lambda 6 = 1 - sum(e^2)/Vx = 1-sum(1-r^2(smc))/Vx.}

The squared multiple correlation is a lower bound for the item communality and as the number of items increases, becomes a better estimate.

G6 is also sensitive to lumpyness in the test and should not be taken as a measure of unifactorial structure.  For lumpy tests, it will be greater than alpha.  For tests with equal item loadings, alpha > G6, but if the loadings are unequal or if there is a general factor, G6 > alpha. alpha is a generalization of an earlier estimate of reliability for tests with dichotomous items developed by Kuder and Richardson, known as KR20, and a shortcut approximation, KR21. (See Revelle, in prep).

Alpha and G6 are both positive functions of the number of items in a test as well as the average intercorrelation of the items in the test.  When calculated from the item variances and total test variance, as is done here, raw alpha is sensitive to differences in the item variances. Standardized alpha is based upon the correlations rather than the covariances.  

A useful index of the quality of the test that is linear with the number of items and the average correlation is the Signal/Noise ratio where \deqn{s/n = \frac{n \bar{r}}{1-\bar{r}}}{s/n = n r/(1-r)}  (Cronbach and Gleser, 1964; Revelle and Condon (in press)).


More complete reliability analyses of a single scale can be done using the \code{\link{omega}} function which finds \eqn{\omega_h}{omega_hierchical} and \eqn{\omega_t}{omega_total} based upon a hierarchical factor analysis.  

Alternative functions \code{\link{score.items}} and  \code{\link{cluster.cor}} will also score multiple scales and report more useful statistics. ``Standardized" alpha is calculated from the inter-item correlations and will differ from raw alpha. 

Four alternative item-whole correlations are reported, three are conventional, one unique.  raw.r is the correlation of the item with the entire scale, not correcting for item overlap. std.r is the correlation of the item with the entire scale, if each item were standardized.  r.drop is the correlation of the item with the scale composed of the remaining items.  Although each of these are conventional statistics, they have the disadvantage that a) item overlap inflates the first and b) the scale is different for each item when an item is dropped. Thus, the fourth alternative, r.cor, corrects for the item overlap by subtracting the item variance but then replaces this with the best estimate of common variance, the smc. This is similar to a suggestion by Cureton (1966).

If some items are to be reversed keyed then they can be specified by either item name or by item location.  (Look at the 3rd and 4th examples.)  Automatic reversal can also be done, and this is based upon the sign of the loadings on the first principal component (Example 5).  This requires the check.keys option to be TRUE.  Previous versions defaulted to have check.keys=TRUE, but some users complained that this made it too easy to find alpha without realizing that some items had been reversed (even though a warning was issued!).  Thus, I have set the default to be check.keys=FALSE with a warning that some items need to be reversed (if this is the case).  To suppress these warnings, set warnings=FALSE.  

Scores are based upon the simple averages (or totals) of the items scored. Thus, if some items are missing, the scores reflect just the items answered.  This is particularly problematic if using total scores (with the cumulative=TRUE option).  To impute missing data using either means or medians, use the \code{\link{scoreItems}} function.   Reversed items are subtracted from the maximum + minimum item response for all the items.

When using raw data, standard errors for the raw alpha are calculated using equation 2 and 3 from Duhhachek and Iacobucci (2004).  This is problematic because some simulations suggest these values are too small.  It is probably better to use bootstrapped values.

\code{\link{alpha.ci}} finds confidence intervals using  the Feldt et al. (1987) procedure.  This procedure does not consider the internal structure of the test the way that the Duhhachek and Iacobucci (2004) procedure does.  That is, the latter considers the variance of the covariances, while the Feldt procedure is based upon just the mean covariance.

Bootstrapped resamples are found if n.iter > 1.  These are returned as the boot object.  They may be plotted or described.  The 2.5\% and 97.5\% values are returned in the boot.ci object.

}
\value{
  \item{total }{a list containing}
  \item{raw_alpha}{alpha based upon the covariances}
  \item{std.alpha}{The standarized alpha based upon the correlations}
  \item{G6(smc)}{Guttman's Lambda 6 reliability}
  \item{average_r}{The average interitem correlation}
  \item{median_r}{The median interitem correlation}
  \item{mean}{For data matrices, the mean of the scale formed by averaging or summing the items (depending upon the cumulative option)}
  \item{sd}{For data matrices, the standard deviation of the total score}
  \item{alpha.drop }{A data frame with all of the above for the case of each item being removed one by one.}
  \item{item.stats}{A data frame including}
  \item{n}{number of complete cases for the item}
  \item{raw.r}{The correlation of each item with the total score, not corrected for item overlap.}
  \item{std.r}{The correlation of each item with the total score (not corrected for item overlap) if the items were all standardized}
  \item{r.cor}{Item whole correlation corrected for item overlap and scale reliability}
  \item{r.drop}{Item whole correlation for this item against the scale without this item}
  \item{mean}{for data matrices, the mean of each item}
  \item{sd}{For data matrices, the standard deviation of each item}
  \item{response.freq}{For data matrices, the frequency of each item response (if less than 20)}
  \item{scores}{Scores are by default simply the average response for all items that a participant took. If cumulative=TRUE, then these are sum scores. Note, this is dangerous if there are lots of missing values.}  
  item{boot.ci}{The lower, median, and upper ranges of the 95\% confidence interval based upon the bootstrap.}
  \item{boot}{a 6 column by n.iter matrix of boot strapped resampled values}
  \item{Unidim}{An index of unidimensionality}
  \item{Fit}{The fit of the off diagonal matrix}
  
}
\references{
Cronbach, L.J. (1951) Coefficient alpha and the internal strucuture of tests.  Psychometrika, 16, 297-334.

Cureton, E. (1966). Corrected item-test correlations. Psychometrika, 31(1):93-96.

Cronbach, L.J. and Gleser G.C.  (1964)The signal/noise ratio in the comparison of reliability coefficients. Educational and Psychological Measurement, 24 (3) 467-480. 

Duhachek, A. and Iacobucci, D. (2004). Alpha's standard error (ase): An accurate and precise confidence interval estimate. Journal of Applied Psychology, 89(5):792-808.

Feldt, L. S., Woodruff, D. J., & Salih, F. A. (1987). Statistical inference for coefficient alpha. Applied Psychological Measurement (11) 93-103.

Guttman, L. (1945). A basis for analyzing test-retest reliability. Psychometrika, 10 (4), 255-282. 

Revelle, W.  (in preparation) An introduction to psychometric theory with applications in {R}. Springer.  (Available online at \url{https://personality-project.org/r/book/}). 

Revelle, W. Hierarchical Cluster Analysis and the Internal Structure of Tests. Multivariate Behavioral Research, 1979, 14, 57-74.

Revelle, W. and Condon, D.M. (2019) Reliability from alpha to omega: A tutorial.  Psychological Assessment, 31, 12, 1395-1411. https://doi.org/10.1037/pas0000754.  \url{https://psyarxiv.com/2y3w9/} Preprint available from PsyArxiv 

Revelle, W. and Condon, D.M. (2018) Reliability.  In Irwing, P., Booth, T. and Hughes, D. (Eds). the Wiley-Blackwell Handbook of Psychometric Testing: A multidisciplinary reference on survey, scale, and test development.

Revelle, W. and Zinbarg, R. E. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma. Psychometrika, 74 (1) 1145-154. 
}
\author{William Revelle }

\note{By default, items that correlate negatively with the overall scale will be reverse coded.  This option may be turned off by setting check.keys = FALSE.  If items are reversed, then each item is subtracted from the  minimum item response + maximum item response where  min and max are taken over all items.  Thus, if the items intentionally differ in range, the scores will be off by a constant.  See \code{\link{scoreItems}} for a solution. 

Two item level statistics are worth comparing: the mean interitem r and the median interitem r.  If these differ very much, that is a sign that the scale is not particularly homogeneous.  

Variables without variance do not contribute to reliability but do contribute to total score. They are dropped with a warning that they had no variance and were thus dropped.  However the scores found still include these values in the calculations.


If the data have been preprocessed by the dplyr package, a strange error can occur.  alpha expects either data.frames or matrix input. data.frames returned by dplyr have had three extra classes added to them which causes alpha to break.  The solution is merely to change the class of the input to "data.frame".  

Two experimental measures of Goodness of Fit are returned in the output: Unidim and Fit. They are not printed or displayed, but are available for analysis. The first is an index of how well the modeled average correlations actually reproduce the original correlation matrix.  The second is how well the modeled correlations reproduce the off diagonal elements of the matrix. Both are indices of squared residuals compared to the squared original correlations.  These two measures are under development and might well be modified or dropped in subsequent versions.}

\seealso{  \code{\link{omega}}, \code{\link{ICLUST}}, \code{\link{guttman}}, \code{\link{scoreItems}}, \code{\link{cluster.cor}}}
\examples{
set.seed(42) #keep the same starting values
#four congeneric measures
r4 <- sim.congeneric()
alpha(r4)
#nine hierarchical measures -- should actually use omega
r9 <- sim.hierarchical()
alpha(r9)

# examples of two independent factors that produce reasonable alphas
#this is a case where alpha is a poor indicator of unidimensionality

two.f <- sim.item(8)
#specify which items to reverse key by name
 alpha(two.f,keys=c("V3","V4","V5","V6"))
 cov.two <- cov(two.f)
 alpha(cov.two,check.keys=TRUE)
 #automatic reversal base upon first component
alpha(two.f,check.keys=TRUE)    #note that the median is much less than the average R
#this suggests (correctly) that the 1 factor model is probably wrong 
#an example with discrete item responses  -- show the frequencies
items <- sim.congeneric(N=500,short=FALSE,low=-2,high=2,
        categorical=TRUE) #500 responses to 4 discrete items with 5 categories
a4 <- alpha(items$observed)  #item response analysis of congeneric measures
a4
#summary just gives Alpha
summary(a4)
}

\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line
