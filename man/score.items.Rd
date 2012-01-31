\name{score.items}
\alias{score.items}
\alias{response.frequencies}
\title{ Score item composite scales and find Cronbach's alpha, Guttman lambda 6 and item whole correlations }
\description{Given a matrix or data.frame of k keys for m items (-1, 0, 1), and a matrix or data.frame of items scores for m items and n people, find the sum scores or average scores for each person and each scale.  In addition, report Cronbach's alpha, Guttman's Lambda 6,  the average r, the scale intercorrelations, and the item by scale correlations (raw and corrected for item overlap).  Replace missing values with the item median or mean if desired. Will adjust scores for reverse scored items.  See \code{\link{make.keys}} for a convenient way to make the keys file.  If the input is a square matrix, then it is assumed that the input is a covariance or correlation matix and scores are not found, but the item statistics are reported. (Similar functionality to \code{\link{cluster.cor}}). \code{\link{response.frequencies}} reports the frequency of item endorsements fore each response category for polytomous or multiple choice items.  
}
\usage{
score.items(keys, items, totals = FALSE, ilabels = NULL,missing=TRUE, impute="median",  min = NULL, max = NULL, digits = 2)
response.frequencies(items,max=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keys}{A matrix or dataframe of -1, 0, or 1 weights for each item on each scale. May be created by hand, or by using \code{\link{make.keys}} }
  \item{items}{ Matrix or dataframe of raw item scores}
  \item{totals}{ if TRUE  find total scores, if FALSE (default), find average scores }
  \item{ilabels}{ a vector of item labels. }
  \item{missing}{missing = TRUE is the normal case and data are imputed according to the impute option.  missing=FALSE, only complete cases are scored.}
  \item{impute}{impute="median" replaces missing values with the item median, impute = "mean" replaces values with the mean response. impute="none" the subject's scores are based upon the average of the keyed, but non missing scores.  }
  \item{min}{May be specified as minimum item score allowed, else will be calculated from data}
  \item{max}{May be specified as maximum item score allowed, else will be calculated from data.  Alternatively, in response frequencies, it is maximum number  of alternative responses to count. }
  \item{digits}{ Number of digits to report }

  
}
\details{The process of finding sum or average scores for a set of scales given a larger set of items is a typical problem in psychometric research.  Although the structure of scales can be determined from the item intercorrelations, to find scale means, variances, and do further analyses, it is typical to find scores based upon the sum or the average item score.  For some strange reason, personality scale scores are typically given as totals, but attitude scores as averages.  The default for score.items is the average as it would seem to make more sense to report scale scores in the metric of the item.  

Various estimates of scale reliability include ``Cronbach's alpha", Guttman's Lambda 6, and the average interitem correlation.  For k = number of items in a scale, and av.r = average correlation between items in the scale, alpha = k * av.r/(1+ (k-1)*av.r).  Thus, alpha is an increasing function of test length as well as the test homeogeneity.  

Surprisingly, 106 years after Spearman (1904) introduced the concept of reliability to psychologists, there are still multiple approaches for measuring it. Although very popular, Cronbach's \eqn{\alpha}  (1951) underestimates the reliability of a test and over estimates the first factor saturation.

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

G6 is also sensitive to lumpyness in the test and should not be taken as a measure of unifactorial structure.  For lumpy tests, it will be greater than alpha.  For tests with equal item loadings, alpha > G6, but if the loadings are unequal or if there is a general factor, G6 > alpha. Although it is normal when scoring just a single scale to calculate G6 from just those items within the scale, logically it is appropriate to estimate an item reliability from all items available.  This is done here and is labeled as G6* to identify the subtle difference.

Alpha and G6* are both positive functions of the number of items in a test as well as the average intercorrelation of the items in the test.  When calculated from the item variances and total test variance, as is done here, raw alpha is sensitive to differences in the item variances. Standardized alpha is based upon the correlations rather than the covariances.  alpha is a generalization of an earlier estimate of reliability for tests with dichotomous items developed by Kuder and Richardson, known as KR20, and a shortcut approximation, KR21. (See Revelle, in prep).

More complete reliability analyses of a single scale can be done using the \code{\link{omega}} function which finds \eqn{\omega_h}{omega_hierchical} and \eqn{\omega_t}{omega_total} based upon a hierarchical factor analysis.  Alternative estimates of the Greatest Lower Bound for the reliability are found in the \code{\link{guttman}} function. 


Alpha is a poor estimate of the general factor saturation of a test (see Revelle and Zinbarg, 2009; Zinbarg et al., 2005) for it can seriously overestimate the size of a general factor, and a better but not perfect estimate of total test reliability because it underestimates total reliability. None the less, it is a useful statistic to report. 

Correlations between scales are attenuated by a lack of reliability.  Correcting correlations for reliability (by dividing by the square roots of the reliabilities of each scale) sometimes help show structure.  

By default, missing values are replaced with the corresponding median value for that item.  Means can be used instead (impute="mean"), or subjects with missing data can just be dropped (missing = FALSE).  For data with a great deal of missingness, yet another option is to just find the average of the available responses (impute="none").  This is useful for findings means for scales for the SAPA project where most scales are estimated from random sub samples of the items from the scale. In this case, the alpha reliabilities are seriously overinflated because they are based upon the total number of items in each scale.  The "alpha observed" values are based upon the average number of items answered in each scale using the standard form for alpha a function of inter-item correlation and number of items.


}
\value{ 
 \item{scores }{Sum or average scores for each subject on the k scales}
  \item{alpha }{Cronbach's coefficient alpha.  A simple (but non-optimal) measure of the internal consistency of a test. See also beta and omega. Set to 1 for scales of length 1. }
  \item{av.r}{The average correlation within a scale, also known as alpha 1, is a useful index of the internal consistency of a domain. Set to 1 for scales with 1 item.}
  \item{G6}{Guttman's Lambda 6 measure of reliability}
  \item{n.items}{Number of items on each scale}
   \item{item.cor}{The correlation of each item with each scale.  Because this is not corrected for item overlap, it will overestimate the amount that an item correlates with the other items in a scale.}
  \item{cor}{The intercorrelation of all the scales}
  \item{corrected}{The correlations of all scales (below the diagonal), alpha on the diagonal, and the unattenuated correlations (above the diagonal)}
  \item{item.corrected}{The item by scale correlations for each item, corrected for item overlap by replacing the item variance with the smc for that item}
 \item{response.freq}{The response frequency (based upon number of non-missing responses) for each alternative.}
 \item{missing}{How many items were not answered  for each scale }
 \item{num.ob.item}{The average number of items with responses on a scale.  Used in calculating the alpha.observed-- relevant for SAPA type data structures.}
}

\note{It is important to recognize in the case of massively missing data (e.g., data from a Synthetic Aperture Personality Assessment study where perhaps only 10-50\% of the items per scale are given to any one subject)) that the number of items per scale, and hence standardized alpha, is not the nominal value and hence alpha will be overestimated.  For this case (impute="none"), an additional alpha (alpha.ob) is reported. 

An alternative to classical test theory scoring is to use \code{\link{score.irt}} to find score estimates based upon Item Response Theory. This is particularly useful in the case of SAPA data which tend to be massively missing.  

} 
\references{
Revelle, W. (in preparation) An introduction to psychometric theory with applications in R.  \url{http://personality-project.org/r/book}


Revelle W. and R.E. Zinbarg. (2009) 
Coefficients alpha, beta, omega and the glb: comments on Sijtsma.
Psychometrika, 74(1):145-154.


Zinbarg, R. E., Revelle, W., Yovel, I. and Li, W. (2005) Cronbach's alpha, Revelle's beta, and McDonald's omega h, Their relations with each other and two alternative conceptualizations of reliability, Psychometrika, 70, 123-133.
	
}
\author{ William Revelle }

\seealso{\code{\link{make.keys}} for a convenient way to create the keys file,  \code{\link{score.multiple.choice}} for multiple choice items, \cr \code{\link{alpha}}, \code{\link{correct.cor}}, \code{\link{cluster.cor}} , \code{\link{cluster.loadings}}, \code{\link{omega}}, \code{\link{guttman}} for item/scale analysis.

In addition, the \code{\link{irt.fa}} function provides an alternative way of examining the structure of a test and emphasizes item response theory approaches to the information returned by each item and the total test. Associated with these IRT parameters is the \code{\link{score.irt}} function for finding IRT based scores.}
\examples{

#see  the example including the bfi data set
data(bfi)
 keys.list <- list(agree=c(-1,2:5),conscientious=c(6:8,-9,-10),extraversion=c(-11,-12,13:15),neuroticism=c(16:20),openness = c(21,-22,23,24,-25))
 keys <- make.keys(28,keys.list,item.labels=colnames(bfi))
 scores <- score.items(keys,bfi,min=1,max=6)
 summary(scores)
 #to get the response frequencies, we need to not use the age variable
 scores <- score.items(keys[1:27,],bfi[1:27],min=1,max=6)
 scores
 
 #compare this output to that for the impute="none" option for SAPA type data
 #first make many of the items missing in a missing pattern way
 missing.bfi <- bfi
 missing.bfi[1:1000,3:8] <- NA
 missing.bfi[1001:2000,c(1:2,9:10)] <- NA
 scores <- score.items(keys,missing.bfi,impute="none",min=1,max=6)
 scores

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}

