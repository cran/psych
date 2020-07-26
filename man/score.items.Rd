\name{scoreItems}
\alias{scoreItems}
\alias{scoreFast}
\alias{scoreVeryFast}
\alias{score.items}
\alias{response.frequencies}
\alias{responseFrequency}
\title{ Score item composite scales and find Cronbach's alpha, Guttman lambda 6 and item whole correlations }
\description{Given a data.frame or matrix of n items and N observations and a list of the direction to score them (a keys.list with k keys)  find the sum scores or average scores for each person and each scale.  In addition, report Cronbach's alpha, Guttman's Lambda 6,  the average r, the scale intercorrelations, and the item by scale correlations (raw and corrected for item overlap).  Replace missing values with the item median or mean if desired. Items may be keyed 1 (score it), -1 ) (reverse score it), or 0 (do not score it). Negatively keyed items will be reverse scored.  Although prior versions used a keys matrix, it is now recommended to just use a list of scoring keys. See \code{\link{make.keys}} for a convenient way to make the keys file.  If the input is a square matrix, then it is assumed that the input is a covariance or correlation matix and scores are not found, but the item statistics are reported. (Similar functionality to \code{\link{cluster.cor}}). \code{\link{response.frequencies}} reports the frequency of item endorsements fore each response category for polytomous or multiple choice items. \code{\link{scoreFast}} and   \code{\link{scoreVeryFast}} just find sums/mean scores and do not report reliabilities.  Much faster for large data sets.
}
\usage{
scoreItems(keys, items, totals = FALSE, ilabels = NULL,missing=TRUE, impute="median",
      delete=TRUE, min = NULL, max = NULL, digits = 2,n.obs=NULL,select=TRUE)
score.items(keys, items, totals = FALSE, ilabels = NULL,missing=TRUE, impute="median",
      delete=TRUE, min = NULL, max = NULL, digits = 2,select=TRUE) 
scoreFast(keys, items, totals = FALSE, ilabels = NULL,missing=TRUE, impute="none",
       delete=TRUE, min = NULL, max = NULL,count.responses=FALSE, digits = 2)
scoreVeryFast(keys,items,totals=FALSE, min=NULL,max=NULL,count.responses=FALSE)
response.frequencies(items,max=10,uniqueitems=NULL)
responseFrequency(items,max=10,uniqueitems=NULL)
}

\arguments{
  \item{keys}{A list of scoring keys or a matrix or dataframe of -1, 0, or 1 weights for each item on each scale which may be created by hand, or by using \code{\link{make.keys}}. Just using a list of scoring keys (see example) is probably more convenient.}
  \item{items}{ Matrix or dataframe of raw item scores}
  \item{totals}{ if TRUE  find total scores, if FALSE (default), find average scores }
  \item{ilabels}{ a vector of item labels. }
  \item{missing}{missing = TRUE is the normal case and data are imputed according to the impute option.  missing=FALSE, only complete cases are scored.}
  \item{impute}{impute="median" replaces missing values with the item medians, impute = "mean" replaces values with the mean response. impute="none" the subject's scores are based upon the average of the keyed, but non missing scores. impute = "none" is probably more appropriate for a large number of missing cases (e.g., SAPA data).  }
  \item{delete}{if delete=TRUE, automatically delete items with no variance (and issue a warning)}
  \item{min}{May be specified as minimum item score allowed, else will be calculated from data.  min and max should be specified if items differ in their possible minima or maxima.  See notes for details.}
  \item{max}{May be specified as maximum item score allowed, else will be calculated from data.  Alternatively, in response frequencies, it is maximum number  of alternative responses to count. }
  \item{uniqueitems}{If specified, the set of possible unique response categories}
  \item{digits}{ Number of digits to report for mean scores }
  \item{n.obs}{If scoring from a correlation matrix, specify the number of subjects allows for the calculation of the confidence intervals for alpha.}
  \item{select}{By default, just find the statistics of those items that are included in scoring keys.  This allows scoring of data sets that have bad data for some items that are not included in the scoring keys. This also speeds up the scoring of small subsets of item from larger data sets.}
\item{count.responses}{If TRUE, report the number of items/scale answered for each subject.}

  
}
\details{The process of finding sum or average scores for a set of scales given a larger set of items is a typical problem in applied psychometrics and in psychometric research.  Although the structure of scales can be determined from the item intercorrelations, to find scale means, variances, and do further analyses, it is typical to find scores based upon the sum or the average item score.  For some strange reason, personality scale scores are typically given as totals, but attitude scores as averages.  The default for scoreItems is the average as it would seem to make more sense to report scale scores in the metric of the item.  

When scoring more than one scale, it is convenient to have a list of the items on each scale and the direction to score the items.  This may be converted to a keys.matrix using \code{\link{make.keys}} or may be entered as a keys.list directly.  

Various estimates of scale reliability include ``Cronbach's alpha", Guttman's Lambda 6, and the average interitem correlation.  For k = number of items in a scale, and av.r = average correlation between items in the scale, alpha = k * av.r/(1+ (k-1)*av.r).  Thus, alpha is an increasing function of test length as well as the test homeogeneity.  

Surprisingly, more than a century after Spearman (1904) introduced the concept of reliability to psychologists, there are still multiple approaches for measuring it. Although very popular, Cronbach's \eqn{\alpha}  (1951) underestimates the reliability of a test and over estimates the first factor saturation.

\eqn{\alpha}{alpha} (Cronbach, 1951) is the same as Guttman's  \eqn{\lambda_3}{lambda3} (Guttman, 1945) and may be found by
\deqn{
\lambda_3 =  \frac{n}{n-1}\Bigl(1 - \frac{tr(\vec{V})_x}{V_x}\Bigr) = \frac{n}{n-1} \frac{V_x - tr(\vec{V}_x)}{V_x} = \alpha
}{Lambda 3 = (n)/(n-1)(1-tr(Vx)/(Vx)  = (n)/(n-1)(Vx-tr(Vx)/Vx  = alpha}

Perhaps because it is so easy to calculate and is available in most commercial programs, alpha is without doubt the most frequently reported measure of internal consistency reliability. Alpha is the mean of all possible spit half reliabilities (corrected for test length).  For a unifactorial test, it is a reasonable estimate of the first factor saturation, although if the test has any microstructure (i.e., if it is ``lumpy") coefficients \eqn{\beta}{beta} (Revelle, 1979; see \code{\link{ICLUST}}) and \eqn{\omega_h}{omega_hierchical} (see \code{\link{omega}})  (McDonald, 1999; Revelle and Zinbarg, 2009) are more appropriate estimates of the general factor saturation.  \eqn{\omega_t}{omega_total} (see \code{\link{omega}}) is a better estimate of the reliability of the total test.  

Guttman's Lambda 6 (G6) considers the amount of variance in each item that can be accounted for the linear regression of all of the other items (the squared multiple correlation or smc), or more precisely, the variance of the errors, \eqn{e_j^2},  and is
\deqn{
\lambda_6 = 1 - \frac{\sum e_j^2}{V_x} = 1 - \frac{\sum(1-r_{smc}^2)}{V_x}
.}{lamada 6 = 1 - sum(e^2)/Vx = 1-sum(1-r^2(smc))/Vx.}

The squared multiple correlation is a lower bound for the item communality and as the number of items increases, becomes a better estimate.

G6 is also sensitive to lumpyness in the test and should not be taken as a measure of unifactorial structure.  For lumpy tests, it will be greater than alpha.  For tests with equal item loadings, alpha > G6, but if the loadings are unequal or if there is a general factor, G6 > alpha. Although it is normal when scoring just a single scale to calculate G6 from just those items within the scale, logically it is appropriate to estimate an item reliability from all items available.  This is done here and is labeled as G6* to identify the subtle difference.

Alpha and G6* are both positive functions of the number of items in a test as well as the average intercorrelation of the items in the test.  When calculated from the item variances and total test variance, as is done here, raw alpha is sensitive to differences in the item variances. Standardized alpha is based upon the correlations rather than the covariances.  alpha is a generalization of an earlier estimate of reliability for tests with dichotomous items developed by Kuder and Richardson, known as KR20, and a shortcut approximation, KR21. (See Revelle, in prep; Revelle and Condon, in press.).

A useful index is the ratio of reliable variance to unreliable variance and is known as the Signal/Noise ratio.  This is just \deqn{s/n = \frac{n \bar{r}}{1-n \bar{r}}}{s/n = n r/(1-nr)}  (Cronbach and Gleser, 1964; Revelle and Condon (in press)).

Standard errors for unstandardized alpha are reported using the formula from Duhachek and Iacobucci (2005).

More complete reliability analyses of a single scale can be done using the \code{\link{omega}} function which finds \eqn{\omega_h}{omega_hierchical} and \eqn{\omega_t}{omega_total} based upon a hierarchical factor analysis.  Alternative estimates of the Greatest Lower Bound for the reliability are found in the \code{\link{guttman}} function. 


Alpha is a poor estimate of the general factor saturation of a test (see Revelle and Zinbarg, 2009; Zinbarg et al., 2005) for it can seriously overestimate the size of a general factor, and a better but not perfect estimate of total test reliability because it underestimates total reliability. None the less, it is a common statistic to report. In general, the use of alpha should be discouraged and the use of more appropriate estimates (\eqn{\omega_h}{omega_hierchical} and \eqn{\omega_t}{omega_total}) should be encouraged. 

Correlations between scales are attenuated by a lack of reliability.  Correcting correlations for reliability (by dividing by the square roots of the reliabilities of each scale) sometimes help show structure.  This is done in the scale intercorrelation matrix with raw correlations below the diagonal and unattenuated correlation above the diagonal.


There are several alternative ways to treat missing values.  By default, missing values are replaced with the corresponding median value for that item.  Means can be used instead (impute="mean"), or subjects with missing data can just be dropped (missing = FALSE).  For data with a great deal of missingness, yet another option is to just find the average of the available responses (impute="none").  This is useful for findings means for scales for the SAPA project (see \url{https://sapa-project.org}) where most scales are estimated from random sub samples of the items from the scale. In this case, the alpha reliabilities are seriously overinflated because they are based upon the total number of items in each scale.  The "alpha observed" values are based upon the average number of items answered in each scale using the standard form for alpha a function of inter-item correlation and number of items. 

Using the impute="none" option as well as asking for totals (totals="TRUE") will be done, although a warning will be issued because scores will now reflect the number of items responded to much more than the actual pattern of responses.

The number of missing responses for each person for each scale is reported in the missing object.  One possibility is to drop scores just for those scales with missing responses.  This may be done adding the code:

scores$scores[scores$missing >0] <- NA 

This is shown in the last example.

Note that the default for scoreItems is to impute missing items with their median, but the default for scoreFAst is to not impute but must return the scale scores based upon the mean or total value for the items scored.




\code{\link{scoreItems}} can be applied to correlation matrices to find just the reliability statistics.  This will be done automatically if the items matrix is symmetric.  

\code{\link{scoreFast}} just finds the scores (with or without imputation) and does not report other statistics.  It is much faster!
\code{\link{scoreVeryFast}} is even more stripped down, no imputation, just scores based upon the observed data. No statistics.
}
\value{ 
 \item{scores }{Sum or average scores for each subject on the k scales}
  \item{alpha }{Cronbach's coefficient alpha.  A simple (but non-optimal) measure of the internal consistency of a test. See also beta and omega. Set to 1 for scales of length 1. }
  \item{av.r}{The average correlation within a scale, also known as alpha 1, is a useful index of the internal consistency of a domain. Set to 1 for scales with 1 item.}
  \item{G6}{Guttman's Lambda 6 measure of reliability}
  \item{G6*}{A generalization of Guttman's Lambda 6 measure of reliability using all the items to find the smc.}
  \item{n.items}{Number of items on each scale}
   \item{item.cor}{The correlation of each item with each scale.  Because this is not corrected for item overlap, it will overestimate the amount that an item correlates with the other items in a scale.}
  \item{cor}{The intercorrelation of all the scales based upon the interitem correlations (see note for why these differ from the correlations of the observed scales themselves).}
  \item{corrected}{The correlations of all scales (below the diagonal), alpha on the diagonal, and the unattenuated correlations (above the diagonal)}
  \item{item.corrected}{The item by scale correlations for each item, corrected for item overlap by replacing the item variance with the smc for that item}
 \item{response.freq}{The response frequency (based upon number of non-missing responses) for each alternative.}
 \item{missing}{How many items were not answered  for each scale }
 \item{num.ob.item}{The average number of items with responses on a scale.  Used in calculating the alpha.observed-- relevant for SAPA type data structures.}
}

\note{It is important to recognize in the case of massively missing data (e.g., data from a Synthetic Aperture Personality Assessment (\url{https://sapa-project.org}) or the International Cognitive Ability Resources (\url{https://icar-project.org})) study where perhaps only 10-50\% of the items per scale are given to any one subject)) that the number of items per scale, and hence standardized alpha, is not the nominal value and hence alpha of the observed scales will be overestimated.  For this case (impute="none"), an additional alpha (alpha.ob) is reported. 

More importantly in this case of massively missing data, there is a difference between the correlations of the composite scales based upon the correlations of the items and the correlations of the scored scales based upon the observed data.  That is, the cor object will have correlations as if all items had been given, while the correlation of the scores object will reflect the actual correlation of the scores.  For SAPA data, it is recommended to use the cor object.  Confidence of these correlations may be found using the \code{\link{cor.ci}} function.

Further note that the inter-scale correlations are based upon the correlations of scales formed from the covariance matrix of the items.  This will differ from the correlation  of scales based upon the correlation of the items.  Thus, although  \code{\link{scoreItems}} will produce reliabilities and intercorrelations from either the raw data or from a correlation matrix, these values will differ slightly.  In addition, with a great deal of missing data, the scale intercorrelations will differ from the correlations of the scores produced, for the latter will be attenuated. 

An alternative to classical test theory scoring is to use \code{\link{scoreIrt}} to find score estimates based upon Item Response Theory. This is particularly useful in the case of SAPA data which tend to be massively missing.  It is also useful  to find scores based upon polytomous items following a factor analysis of the polychoric correlation matrix (see \code{\link{irt.fa}}).  However, remember that  this only makes sense if the items are unidimensional. That is to say, if forming item composites from (e.g., \code{\link{bestScales}}), that are empirically derived, they will necessarily have a clear factor structure and the IRT based scoring does not make sense.

When reverse scoring items from a set where items differ in their possible minima or maxima, it is important to specify the min and max values.  Items are reversed by subtracting them from max + min.  Thus, if items range from 1 to 6, items are reversed by subtracting them from 7.  But, if the data set includes other variables, (say an id field) that far exceeds the item min or max, then the max id will incorrectly be used to reverse key.  min and max can either be single values, or vectors for all items. Compare two examples of scoring the \code{\link[psychTools]{bfi}} data set.

If scales are formed with overlapping items, then the correlations of the scales will be seriously inflated.  \code{\link{scoreOverlap}} will adjust the correlations for this overlap.

Yet another possibility for scoring large data sets is to ignore all the reliability calculations and just find the scores.  This may be done using \code{\link{scoreFast}} or 
\code{\link{scoreVeryFast}}.  These two functions just find mean scores (\code{\link{scoreVeryFast}} without imputation) or will do imputation if desired \code{\link{scoreFast}}. For 200K cases on 1000 variables with 11 scales, \code{\link{scoreVeryFast}} took 4.7 seconds on a Mac PowerBook with a 2.8GHZ Intel I7 and \code{\link{scoreFast}} took 23.2 seconds.  \code{\link{scoreIrt.1pl}} for the same problem took xxx with options("mc.cores"=1) (not parallel processing) and 1259 seconds with  options("mc.cores"=NULL) (implying 2 cores) and with four cores was very slow (probably too much parallel processing). 

Yet one more possibility is to find scores based upon a matrix of weights (e.g. zero order correlations, beta weights, or factor weights.) In this case, scores are simply the product of the weights times the (standardized) items.   If using coefficients from a regression analysis (lm), a column of 1's is added to the data and labeled "(Intercept)" and this is used in the calculation. This is done by \code{\link{scoreWtd}}.
} 
\references{

Cronbach, L.J. and Gleser G.C.  (1964)The signal/noise ratio in the comparison of reliability coefficients. Educational and Psychological Measurement, 24 (3) 467-480. 


Duhachek, A. and Iacobucci, D. (2004). Alpha's standard error (ase): An accurate and precise confidence interval estimate. Journal of Applied Psychology, 89(5):792-808.


McDonald, R. P. (1999). Test theory: A unified treatment. L. Erlbaum Associates, Mahwah, N.J.

Revelle, W. (in preparation) An introduction to psychometric theory with applications in R.  \url{https://personality-project.org/r/book}

Revelle, W. and Condon, D.C. Reliability.  In Irwing, P., Booth, T. and Hughes, D. (Eds). the Wiley-Blackwell Handbook of Psychometric Testing (in press).

Revelle W. and R.E. Zinbarg. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma. Psychometrika, 74(1):145-154.

Zinbarg, R. E., Revelle, W., Yovel, I. and Li, W. (2005) Cronbach's alpha, Revelle's beta, and McDonald's omega h, Their relations with each other and two alternative conceptualizations of reliability, Psychometrika, 70, 123-133.



	
}
\author{ William Revelle }

\seealso{\code{\link{make.keys}} for a convenient way to create the keys file,  \code{\link{score.multiple.choice}} for multiple choice items, \cr \code{\link{alpha}}, \code{\link{correct.cor}}, \code{\link{cluster.cor}} , \code{\link{cluster.loadings}}, \code{\link{omega}}, \code{\link{guttman}} for item/scale analysis.

If scales are formed from overlapping sets of items, their correlations will be inflated.  This is corrected for when using the \code{\link{scoreOverlap}} function which, although it will not produce scores, will report scale intercorrelations corrected for item overlap. 

In addition, the \code{\link{irt.fa}} function provides an alternative way of examining the structure of a test and emphasizes item response theory approaches to the information returned by each item and the total test. Associated with these IRT parameters is the \code{\link{scoreIrt}} function for finding IRT based scores as well as \code{\link{irt.responses}} to show response curves for the alternatives in a multiple choice test.

\code{\link{scoreIrt}} will find both  IRT based estimates as well as average item response scores.  These latter correlate perfectly with those found by scoreItems.  If using a keys matrix, the score.irt results are based upon the item difficulties with the assumption that all items are equally discriminating (effectively a Rasch model).  These scores are probably most useful  in the case of massively missing data because they can take into account the item difficulties. 

\code{\link{scoreIrt.1pl}} finds the item difficulty parameters and then applies a 1 parameter (Rasch like) model.  It chooses items based upon a keys.list.  

}
\examples{

#see  the example including the bfi data set
data(psychTools::bfi)
keys.list <- list(agree=c("-A1","A2","A3","A4","A5"),
  conscientious=c("C1","C2","C3","-C4","-C5"),extraversion=c("-E1","-E2","E3","E4","E5"),
  neuroticism=c("N1","N2","N3","N4","N5"), openness = c("O1","-O2","O3","O4","-O5")) 
  keys <- make.keys(psychTools::bfi,keys.list)  #no longer necessary
 scores <- scoreItems(keys,psychTools::bfi,min=1,max=6)  #using a keys matrix 
 scores <- scoreItems(keys.list,psychTools::bfi,min=1,max=6)  # or just use the keys.list
 summary(scores)
 #to get the response frequencies, we need to not use the age variable
 scores <- scoreItems(keys[1:25,],psychTools::bfi[1:25]) #we do not need to specify min or 
 #max if  there are no values (such as age) outside the normal item range.
 scores
 #The scores themselves are available in the scores$scores object.  I.e.,
 describe(scores$scores)
 
 
 #compare this output to that for the impute="none" option for SAPA type data
 #first make many of the items missing in a missing pattern way
 missing.bfi <- psychTools::bfi
 missing.bfi[1:1000,3:8] <- NA
 missing.bfi[1001:2000,c(1:2,9:10)] <- NA
 scores <- scoreItems(keys.list,missing.bfi,impute="none",min=1,max=6)
 scores
 describe(scores$scores)  #the actual scores themselves
 
 #If we want to delete scales scores for people who did not answer some items for one 
 #(or more) scales, we can do the following:
 
  scores <- scoreItems(keys.list,missing.bfi,totals=TRUE,min=1,max=6) #find total scores
  describe(scores$scores) #note that missing data were replaced with median for the item
  scores$scores[scores$missing > 0] <- NA  #get rid of cases with missing data
  describe(scores$scores)


}
\keyword{ multivariate }
\keyword{models}

