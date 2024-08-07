\name{statsBy}
\alias{statsBy}
\alias{statsBy.boot}
\alias{statsBy.boot.summary}
\alias{faBy}
\title{Find statistics (including correlations) within and between groups for basic multilevel analyses}
\description{When examining data at two levels (e.g., the individual and by some set of grouping variables), it is useful to find basic descriptive statistics (means, sds, ns per group, within group correlations) as well as between group statistics (over all descriptive statistics, and overall between group correlations). Of particular use is the ability to decompose a matrix of correlations at the individual level into correlations within group and correlations between groups. 
}
\usage{
statsBy(data, group, cors = FALSE, cor="cor", method="pearson", use="pairwise", 
	poly=FALSE, na.rm=TRUE,alpha=.05,minlength=5,
	weights=NULL,mean.weights=NULL, min.n=1)
statsBy.boot(data,group,ntrials=10,cors=FALSE,replace=TRUE,method="pearson")
statsBy.boot.summary(res.list,var="ICC2")
faBy(stats, nfactors = 1, rotate = "oblimin", fm = "minres", free = TRUE, all=FALSE,
   min.n = 12,quant=.1, ...)
}
\arguments{
  \item{data}{A matrix or dataframe with rows for subjects, columns for variables.  
  	One of these columns should be the values of a grouping variable.}
  \item{group}{The names or numbers of the variable in data to use as the grouping variables.}
  \item{cors}{Should the results include the correlation matrix within each group?  Default is FALSE.}
  \item{cor}{Type of correlation/covariance to find within groups and between groups.  
  	The default is Pearson correlation.  To find within and between covariances, set cor="cov". 
  	 Although polychoric, tetrachoric, and mixed correlations can be found within groups, 
  	 this does not make sense for the between groups or the pooled within groups.  
  	 In this case, correlations for each group will be as specified, but the between groups 
  	 and pooled within will be Pearson. See the discussion below.}
\item{method}{What kind of correlations should be found (default is Pearson product moment)}
\item{use}{How to treat missing data.  use="pairwise" is the default}
\item{poly}{Find polychoric.tetrachoric correlations within groups if requested.}
\item{na.rm}{Should missing values be deleted (na.rm=TRUE) or should we assume the data are clean?}
\item{alpha}{The alpha level for the confidence intervals for the ICC1 and ICC2, and rwg, rbg}
\item{minlength}{The minimum length to use when abbreviating the labels for confidence intervals}
\item{weights}{If specified, weight the groups by weights when finding the pooled within group correlation.
	  Otherwise weight by sample size.}
\item{mean.weights}{If not NULL, what weights should be use we wait the withingroup variables by some weights vector?}
\item{ntrials}{The number of trials to run when bootstrapping statistics}
\item{replace}{Should the bootstrap be done by permuting the data (replace=FALSE) or sampling with replacement (replace=TRUE)}
\item{res.list}{The results from statsBy.boot may be summarized using boot.stats}
\item{var}{Name of the variable to be summarized from statsBy.boot}
\item{stats}{The output of statsBy}
\item{nfactors}{The number of factors to extract in each subgroup}
\item{rotate}{The factor rotation/transformation}
\item{fm}{The factor method (see \code{\link{fa}} for details)}
\item{free}{Allow the factor solution to be freely estimated for each individual (see note).}
\item{all}{Report individual factor analyses for each group as well as the summary table}
\item{min.n}{The minimum number of within subject cases before we factor analyze it or find group statistics. Only works for one grouping variable.}
\item{quant}{Show the upper and lower quant quantile of the factor loadings in faBy}
\item{...}{Other parameters to pass to the fa function}
}
\details{Multilevel data are endemic in psychological research. In multilevel data, observations are taken on subjects who are nested within some higher level grouping variable.  The data might be experimental (participants are nested within experimental conditions) or observational (students are nested within classrooms, students are nested within college majors.) To analyze this type of data, one uses random effects models or mixed effect models, or more generally, multilevel models.  There are at least two very powerful packages (nlme and multilevel) which allow for complex analysis of hierarchical (multilevel) data structures.  \code{\link{statsBy}} is a much simpler function to give some of the basic descriptive statistics for two level models.  It is meant to supplement true multilevel modeling.

For a group variable (group) for a data.frame or matrix (data), basic descriptive statistics (mean, sd, n) as well as within group correlations (cors=TRUE) are found for each group.  

The amount of variance associated with the grouping variable compared to the total variance is the type 1 IntraClass Correlation (ICC1):
\eqn{ICC1 = (MSb-MSw)/(MSb + MSw*(npr-1))}
where npr is the average number of cases within each group. 

The reliability of the group differences may be found by the ICC2 which reflects how different the means are with respect to the within group variability.  
\eqn{ICC2 = (MSb-MSw)/MSb}.
Because the mean square between is sensitive to sample size, this estimate will also reflect sample size.

Perhaps the most useful part of \code{\link{statsBy}} is that it decomposes the observed correlations between variables into two parts: the within group and the between group correlation. This follows the decomposition of an observed correlation into the pooled correlation within groups (rwg) and the weighted correlation of the means between groups  discussed by Pedazur (1997) and by Bliese in the multilevel package.  

\eqn{r_{xy} = eta_{x_{wg}} * eta_{y_{wg}} * r_{xy_{wg}}  +  eta_{x_{bg}} * eta_{y_{bg}} * r_{xy_{bg}}  }

where \eqn{r_{xy}} is the normal correlation which may be decomposed into a within group and between group correlations \eqn{r_{xy_{wg}}} and \eqn{r_{xy_{bg}}} and eta is the correlation of the data with the within group values, or the group means.

It is important to realize that the within group and between group correlations are independent of each other.  That is to say, inferring from the 'ecological correlation' (between groups) to the lower level (within group) correlation is inappropriate.  However, these between group correlations are still very meaningful, if inferences are made at the higher level.  

There are actually two ways of finding the within group correlations pooled across groups.  We can find the correlations within every group, weight these by the sample size and then report this pooled value (pooled).  This is found if the cors option is set to TRUE.  It is logically  equivalent to doing a sample size weighted meta-analytic correlation.  The other way, rwg, considers the covariances, variances, and thus correlations when each subject's scores are given as deviation score from the group mean.  

If finding tetrachoric, polychoric, or mixed correlations, these two estimates will differ, for the pooled value is the weighted polychoric correlation, but the rwg is the Pearson correlation. 

If the weights parameter is specified, the pooled correlations are found by weighting the groups by the specified weight, rather than sample size.



Confidence values and significance  of  \eqn{r_{xy_{wg}}}, pwg, reflect the pooled number of cases within groups, while  \eqn{r_{xy_{bg}} }, pbg, the number of groups. These are not corrected for multiple comparisons.

\code{\link{withinBetween}} is an example data set of the mixture of within and between group correlations. \code{\link{sim.multilevel}} will generate simulated data with a multilevel structure.

The \code{\link{statsBy.boot}} function will randomize the grouping variable ntrials times and find the statsBy output.  This can take a long time and will produce a great deal of output.  This output can then be summarized for relevant variables using the \code{\link{statsBy.boot.summary}} function specifying the variable of interest.  These two functions are useful in order to find if the mere act of grouping leads to large between group correlations.

Consider the case of the relationship between various tests of ability when the data are grouped by level of education (statsBy(sat.act,"education")) or when affect data are analyzed within and between an affect manipulation (statsBy(flat,group="Film") ). Note in this latter example, that because subjects were randomly assigned to Film condition for the pretest, that the pretest ICC1s cluster around 0. 

\code{\link{faBy}} uses the output of \code{\link{statsBy}} to perform a factor analysis on the correlation matrix within each group. If the free parameter is FALSE, then each solution is rotated towards the group solution (as much as possible).  The output is a list of each factor solution, as well as a summary matrix of loadings and interfactor correlations for  all groups.

}
\value{
\item{means}{The means for each group for each variable.  }
\item{sd}{The standard deviations for each group for each variable.}
\item{n}{The number of cases for each group and for each variable.}
\item{ICC1}{The intraclass correlation reflects the amount of total variance associated with the grouping variable.}
\item{ICC2}{The intraclass correlation (2) reflecting how much the groups means differ.}
\item{ci1}{The confidence intervals for the ICC1}
\item{ci2}{The confidence intervals for the ICC2}
\item{F}{The F from a one-way anova of group means.}
\item{rwg}{The pooled within group correlations.}
\item{ci.wg}{The confidence intervals of the pooled within group correlations.}
\item{rbg}{The sample size weighted between group correlations. }
\item{c.bg}{The confidence intervals of the rbg values}
\item{etawg}{The correlation of the data with the within group values.}
\item{etabg}{The correlation of the data with the group means.}
\item{pbg}{The probability of the between group correlation}
\item{pwg}{The probability of the within group correlation}
\item{r}{In the case that we want the correlations in each group, r is a list of the within group correlations for every group. Set cors=TRUE}
\item{within}{is just another way of displaying these correlations.  within is a matrix which reports the lower off diagonal correlations as one row for each group.}
\item{pooled}{The sample size weighted correlations.  This is just within weighted by the sample sizes. The cors option must be set to TRUE to get this. See the note. }
\item{slope}{The within group slopes --  if using cor="cov" these will be the b weights}

}
\references{ Pedhazur, E.J. (1997) Multiple regression in behavioral research: explanation and prediction.  Harcourt Brace.


}
\author{William Revelle
}
\note{If finding polychoric correlations, the two estimates of the pooled within group correlations will differ, for the pooled value is the weighted polychoric correlation, but the rwg is the Pearson correlation. As of March, 2020 the average is based upon the weighted FisherZ correlation (back transformed) rather than the average correlation. 

The value of rbg (the between group correlation) is the group size weighted correlation.  This is not the same as just finding the correlation of the group means (i.e. cor(means)).


The statsBy.boot function will sometimes fail if sampling with replacement because if the group sizes differ drastically, some groups will be empty.  In this case, sample without replacement.  

The statsBy.boot function can take a long time.  (As I am writing this, I am running 1000 replications of a problem with 64,000 cases and 84 groups.  It is taking about 3 seconds per replication on a MacBook Pro.)

The \code{\link{faBy}} function takes the output of statsBy (with the cors=TRUE option) and then factors each individual subject.  By default, the solutions are organized so that the factors "match" the group solution in terms of their order.  It is also possible to attempt to force the solutions to match by order and also by using the TargetQ rotation function.  (free=FALSE)
}

\seealso{\code{\link{describeBy}} and the functions within the multilevel package.

}
\examples{
#Taken from Pedhazur, 1997
pedhazur <- structure(list(Group = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 
2L), X = c(5L, 2L, 4L, 6L, 3L, 8L, 5L, 7L, 9L, 6L), Y = 1:10), .Names = c("Group", 
"X", "Y"), class = "data.frame", row.names = c(NA, -10L))
pedhazur
ped.stats <- statsBy(pedhazur,"Group")
ped.stats


#Now do this for the sat.act data set
sat.stats <- statsBy(sat.act,c("education","gender"),cors=TRUE)   #group by two grouping variables
print(sat.stats,short=FALSE)
lowerMat(sat.stats$pbg)  #get the probability values

#show means by groups
round(sat.stats$mean)

#Do separate factor analyses for each group
sb.bfi <- statsBy(bfi[1:10], group=bfi$gender,cors=TRUE)
faBy(sb.bfi,1)  #one factor per group 
faBy(sb.bfi,2) #two factors per group



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{ models }
