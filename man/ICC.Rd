\name{ICC}
\alias{ICC}

\title{ Intraclass Correlations (ICC1, ICC2, ICC3 from Shrout and Fleiss) }
\description{The Intraclass correlation is used as a measure of association when studying the reliability of raters.  Shrout and Fleiss (1979) outline 6 different estimates, that depend upon the particular experimental design. All are implemented and given confidence limits.  Uses either aov or lmer depending upon options.  lmer allows for missing values.
}
\usage{
ICC(x,missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE)
}

\arguments{
  \item{x}{a matrix or dataframe of ratings}
  \item{missing}{if TRUE, remove missing data -- work on complete cases only (aov only)}
  \item{alpha}{The alpha level for significance for finding the confidence intervals}
  \item{lmer}{Should we use the lmer function from lme4?  This handles missing data and gives variance components as well. TRUE by default.}
  \item{check.keys}{If TRUE reverse those items that do not correlate with total score. This is not done by default.}
  }
\details{Shrout and Fleiss (1979) consider six cases of reliability of ratings done by k raters on n targets. McGraw and Wong (1996) consider 10, 6 of which are identical to Shrout and Fleiss and  4 are conceptually different but use the same equations as the 6 in Shrout and Fleiss. 

The intraclass correlation is used if raters are all of the same ``class".  That is, there is no logical way of distinguishing them.  Examples include correlations between  pairs of twins, correlations between raters.  If the variables are logically distinguishable (e.g., different items on a test), then the more typical coefficient is based upon the inter-class correlation (e.g., a Pearson r) and a statistic such as \code{\link{alpha}} or \code{\link{omega}} might be used.  alpha and ICC3k are identical.


Where the data are laid out in terms of Rows (subjects) and Columns (rater or tests), the various ICCs are found by the ratio of various estimates of variance components.  In all cases, subjects are taken as varying at random, and the residual variance is also random.  The distinction between models 2 and 3 is whether the judges (items/tests) are seen as random or fixed.  A further distinction is whether the emphasis is upon absolute agreement of the judges, or merely consistency.  

As discussed by Liljequist et al. (2019), McGraw and Wong lay out 5 models which use just three forms of the ICC equations.

Model 1  is a one way model with 

ICC1: Each  target is rated by a different  judge and the judges are selected at random.  
\deqn{
ICC(1,1) = \rho_{1,1} = \frac{\sigma^2_r}{\sigma^2_r + \sigma^2_w} }
(This is a one-way ANOVA fixed effects model and is found by  (MSB- MSW)/(MSB+ (nr-1)*MSW))   

ICC2: A random sample of k judges rate each target.  The measure is one of absolute agreement in the ratings.  
\deqn{
ICC(2,1) = \rho_{2,1} = \frac{\sigma^2_r}{\sigma^2_r + \sigma^2_c +\sigma^2_{rc}  + \sigma^2_e} }

 Found as (MSB- MSE)/(MSB + (nr-1)*MSE + nr*(MSJ-MSE)/nc) 

ICC3: A fixed set of k judges rate each target. There is no generalization to a larger population of judges. 
\deqn{
ICC(3,1) = \rho_{3,1} = \frac{\sigma^2_r}{\sigma^2_r + \sigma^2_c + \sigma^2_e} }


 (MSB - MSE)/(MSB+ (nr-1)*MSE)

Then, for each of these cases, is reliability to be estimated for a single rating or for the average of k ratings?  (The 1 rating case is equivalent to the average intercorrelation, the k rating case to the Spearman Brown adjusted reliability.)

ICC1 is sensitive to differences in means between raters and is a measure of absolute agreement.

ICC2 and ICC3 remove mean differences between judges, but are sensitive to interactions of raters by judges.  The difference between ICC2 and ICC3 is whether raters are seen as fixed or random effects.

ICC1k, ICC2k, ICC3K reflect the means of k raters.   


If using the lmer option, then missing data are allowed.  In addition the lme object returns the variance decomposition.  (This is simliar to  \code{\link{testRetest}} which works on the items from two occasions.

The check.keys option by default reverses items that are negatively correlated with total score.  A message is issued.
}
\value{
  \item{results}{A matrix of 6 rows and 8 columns, including the ICCs, F test, p values, and confidence limits}
  \item{summary}{The anova summary table or the lmer summary table}
  \item{stats}{The anova statistics (converted from lmer if using lmer)}
  \item{MSW}{Mean Square Within based upon the anova}
  \item{lme}{The variance decomposition if using the lmer option}
}

\references{ 
Shrout, Patrick E. and Fleiss, Joseph L. Intraclass correlations: uses in assessing rater reliability.  Psychological Bulletin, 1979, 86, 420-3428.

McGraw, Kenneth O. and Wong, S. P. (1996), Forming inferences about some intraclass correlation coefficients.  Psychological Methods, 1, 30-46. + errata on page 390.

Liljequist David, Elfving Britt and Skavberg Kirsti 
 (2019) Intraclass correlation-A discussion and demonstration of basic features. PLoS ONE 14(7): e0219854. https://doi.org/10.1371/journal.
 
Revelle, W. (in prep) An introduction to psychometric theory with applications in R. Springer. (working draft available at  \url{https://personality-project.org/r/book/}

}
\author{William Revelle }

\note{The results for the Lower and Upper Bounds for ICC(2,k) do not match those of SPSS 9 or 10, but do match the definitions of Shrout and Fleiss.  SPSS  seems to have been using the  formula in McGraw and Wong, but not the errata on p 390.  They seem to have fixed it in more recent releases (15). 

Starting with psych  1.4.2, the confidence intervals are based upon (1-alpha)\% at both tails of the confidence interval.  This is in agreement with Shrout and Fleiss.  Prior to 1.4.2 the confidence intervals were (1-alpha/2)\%. 
However, at some point, this error slipped back again.  It has been fixed in version 1.9.5 (5/21/19).

However, following some very useful comments from Mijke Rhumtulla, I should divide by 2 after all.  This has been fixed (again) 3/05/22 for version 2.2.2.



An occasionally asked question: How do these 6 measures compare to the 10 ICCsreported by SPSS?  The answer is that the 6 Shrout and Fleiss are the same as the 10 by McGraw and Wong.


A very helpful discussion of the formulae used in the calculations is by Liljequist et al, 2019. 


Although M & W distinguish between two way random effects and two way mixed effects for consistency of a single rater, examining their formula show that both of these are S& F ICC(3,1).

Similarly, the S&F ICC(2,1) is called both Two way random effects, absolute agreement and Two way mixed effects, absolute agreement by M&W.

The same confusion occurs with multiple raters (3,k) are both two way random effects and mixed effects for consistency  and (2,k) are both two way random and two way mixed absolute agreement.

As M&W say "The importance of the random-fixed effects distinction is in its effect on the interpretation, but not calculation, of an ICC. Namely, when levels of the column factor are randomly sampled, one can generalize beyond one's data, but not when they are fixed. In either case, however, the value of the ICC is the same," (p 37) 

 
To quote from IBM support:

``Problem

I'm using the SPSS RELIABILITY procedure to compute intraclass correlation coefficients (ICCs), and am trying to make sure I understand the correspondence between the measures available in SPSS and those discussed in Shrout and Fleiss. How do the two sets of models relate to each other?

Resolving The Problem
Shrout and Fleiss (1979) discussed six ICC measures, which consist of pairs of measures for the reliability of a single rating and that of the average of k ratings (where k is the number of raters) for three different models: the one-way random model (called Case 1), the two-way random model (Case 2), and the two-way mixed model (Case 3). The measures implemented in SPSS were taken from McGraw and Wong (1996), who discussed these six measures, plus four others. The additional ones in McGraw and Wong are actually numerically identical to the four measures for the two-way cases discussed by Shrout and Fleiss, differing only in their interpretation.

The correspondence between the measures available in SPSS and those discussed in Shrout and Fleiss is as follows:

S&F(1,1) = SPSS One-Way Random Model, Absolute Agreement, Single Measure

S&F(1,k) = SPSS One-Way Random Model, Absolute Agreement, Average Measure

S&F(2,1) = SPSS Two-Way Random Model, Absolute Agreement, Single Measure

S&F(2,k) = SPSS Two-Way Random Model, Absolute Agreement, Average Measure

S&F(3,1) = SPSS Two-Way Mixed Model, Consistency, Single Measure

S&F(3,k) = SPSS Two-Way Mixed Model, Consistency, Average Measure

SPSS also offers consistency measures for the two-way random case (numerically equivalent to the consistency measures for the two-way mixed case, but differing in interpretation), and absolute agreement measures for the two-way mixed case (numerically equivalent to the absolute agreement measures for the two-way random case, again differing in interpretation)
 
 
 }

\examples{
sf <- matrix(c(
9,    2,   5,    8,
6,    1,   3,    2,
8,    4,   6,    8,
7,    1,   2,    6,
10,   5,   6,    9,
6,   2,   4,    7),ncol=4,byrow=TRUE)
colnames(sf) <- paste("J",1:4,sep="")
rownames(sf) <- paste("S",1:6,sep="")
sf  #example from Shrout and Fleiss (1979)
ICC(sf,lmer=FALSE)  #just use the aov procedure
   
#data(sai)
sai <- psychTools::sai
sai.xray <- subset(sai,(sai$study=="XRAY") & (sai$time==1))
xray.icc <- ICC(sai.xray[-c(1:3)],lmer=TRUE,check.keys=TRUE)
xray.icc
xray.icc$lme  #show the variance components as well
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

