\name{ICC}
\alias{ICC}

\title{ Intraclass Correlations (ICC1, ICC2, ICC3 from Shrout and Fleiss) }
\description{The Intraclass correlation is used as a measure of association when studying the reliability of raters.  Shrout and Fleiss (1979) outline 6 different estimates, that depend upon the particular experimental design. All are implemented and given confidence limits.
}
\usage{
ICC(x,digits=2,alpha=.05)
}

\arguments{
  \item{x}{a matrix or dataframe of ratings}
  \item{digits}{Round the output to digits}
  \item{alpha}{The alpha level for significance for finding the confidence intervals}
}
\details{Shrout and Fleiss (1979) consider six cases of reliability of ratings done by k raters on n targets.

ICC1: Each  target is rated by a different  judge and the judges are selected at random.  (This is a one-way ANOVA fixed effects model and is found by  (MSB- MSW)/(MSB+ (nr-1)*MSW))  

ICC2: A random sample of k judges rate each target.  The measure is one of absolute agreement in the ratings. Found as (MSB- MSE)/(MSB + (nr-1)*MSE + nr*(MSJ-MSE)/nc) 

ICC3: A fixed set of k judges rate each target. There is no generalization to a larger population of judges.  (MSB - MSE)/(MSB+ (nr-1)*MSE)

Then, for each of these cases, is reliability to be estimated for a single rating or for the average of k ratings?  (The 1 rating case is equivalent to the average intercorrelation, the k rating case to the Spearman Brown adjusted reliability.)

ICC1 is sensitive to differences in means between raters and is a measure of absolute agreement.

ICC2 and ICC3 remove mean differences between judges, but are sensitive to interactions of raters by judges.  The difference between ICC2 and ICC3 is whether raters are seen as fixed or random effects.

ICC1k, ICC2k, ICC3K reflect the means of k raters.  

The intraclass correlation is used if raters are all of the same ``class".  That is, there is no logical way of distinguishing them.  Examples include correlations between  pairs of twins, correlations between raters.  If the variables are logically distinguishable (e.g., different items on a test), then the more typical coefficient is based upon the inter-class correlation (e.g., a Pearson r) and a statistic such as \code{\link{alpha}} or \code{\link{omega}} might be used.
}
\value{
  \item{results}{A matrix of 6 rows and 8 columns, including the ICCs, F test, p values, and confidence limits}
  \item{summary}{The anova summary table}
  \item{stats}{The anova statistics}
  \item{MSW}{Mean Square Within based upon the anova}
}
\references{ Shrout, Patrick E. and Fleiss, Joseph L. Intraclass correlations: uses in assessing rater reliability.  Psychological Bulletin, 1979, 86, 420-3428.

McGraw, Kenneth O. and Wong, S. P. (1996), Forming inferences about some intraclass correlation coefficients.  Psychological Methods, 1, 30-46. + errata on page 390.

Revelle, W. (in prep) An introduction to psychometric theory with applications in R. Springer. (working draft available at  \url{http://personality-project.org/r/book/}}

\author{William Revelle }
\note{The results for the Lower and Upper Bounds for ICC(2,k) do not match those of SPSS 9 or 10, but do match the definitions of Shrout and Fleiss.  SPSS  seems to have been using the  formula in McGraw and Wong, but not the errata on p 390.  They seem to have fixed it in more recent releases (15). }

\examples{
sf <- matrix(c(9,    2,   5,    8,
6,    1,   3,    2,
8,    4,   6,    8,
7,    1,   2,    6,
10,   5,   6,    9,
6,   2,   4,    7),ncol=4,byrow=TRUE)
colnames(sf) <- paste("J",1:4,sep="")
rownames(sf) <- paste("S",1:6,sep="")
sf  #example from Shrout and Fleiss (1979)
ICC(sf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

