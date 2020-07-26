\name{testRetest}
\alias{testRetest}
\alias{testReliability}

\title{Find various test-retest statistics, including test, person and item reliability}

\description{Given two presentations of a test, it is straightforward to find the test-retest reliablity, as well as the item reliability and person stability across items.  Using the multi-level structure of the data, it is also possible to do a variance deomposition to find variance components for people, items, time, people x time, people x items, and items x time as well as the residual variance. This leads to various generalizability cofficients.
}
\usage{
testRetest(t1,t2=NULL,keys=NULL,id="id", time=  "time", select=NULL, 
check.keys=TRUE, warnings=TRUE,lmer=TRUE)
}

\arguments{
  \item{t1}{a data.frame or matrix for the first time of measurement.}
  \item{t2}{a data.frame or matrix for the second time of measurement. May be NULL if time is specifed in t1}
  \item{keys}{item names (or locations) to analyze, preface by "-" to reverse score. }
  \item{id}{subject identification codes to match across time}
  \item{time}{The name of the time variable identifying time 1 or 2 if just one data set is supplied. }
  \item{select}{A subset of items to analyze}
  \item{check.keys}{If TRUE will automatically reverse items based upon their correlation with the first principal component.  Will throw a warning when doing so, but some people seem to miss this kind of message.}
  \item{warnings}{If TRUE, then warn when items are reverse scored}
  \item{lmer}{If TRUE, include the lmer variance decomposition. By default, this is true, but this can lead to long times for large data sets. }
}
\details{There are many ways of measuring reliability. Test - Retest is one way.  If the time interval is very short (or immediate), this is known as a dependability correlation, if the time interval is longer, a stability coefficient.  In all cases, this is a correlation between two measures at different time points.  Given the multi-level nature of these data, it is possible to find variance components associated with individuals, time, item, and time by item, etc.  This leads to several different estimates of reliability (see \code{\link{multilevel.reliability}} for a discussion and references).

It is also possible to find the subject reliability across time (this is the correlation across the items at time 1 with time 2 for each subject).  This is a sign of subject reliability (Wood et al, 2017).  Items can show differing amounts of test-retest reliability over time.  Unfortunately, the within person correlation has problems if people do not differ very much across items.  If all items are in the same keyed direction, and measuring the same construct, then the response profile for an individual is essentially flat. This implies that the even with almost perfect reproducibility, that the correlation can actually be negative.  The within person distance (d2) across items is just the mean of the squared differences for each item.  Although highly negatively correlated with the rqq score, this does distinguish between random responders (high dqq and low rqq) from consistent responders with lower variance (low dqq and low rqq). 

Several individual statistics are reported in the scores object.  These can be displayed by using \code{\link{pairs.panels}} for a graphic display of the relationship and ranges of the various measures.

Although meant to decompose the variance for tests with items nested within tests, if just given two tests, the variance components for people and for time will also be shown.  The resulting variance ratio of people to total variance is the intraclass correlation between the two tests.  See also \code{\link{ICC}} for the more general case.
}
\value{
\item{r12}{The time 1 time 2 correlation of scaled scores across time}
\item{alpha}{Guttman's lambda 3 (aka alpha) and lambda 6* (item reliabilities based upon smcs) are found for the scales at times 1 and 2.}
\item{rqq}{The within subject test retest reliability of response patterns over items}
\item{item.stats}{Item reliabilities, item loadings at time 1 and 2, item means at time 1 and time 2}
\item{scores}{A data frame of principal component scores at time 1 and time 2, raw scores from time 1 and time 2, the within person standard deviation for time 1 and time 2, and the rqq  and dqq scores for each subject.  }
\item{xy.df}{If given separate t1 and t2 data.frames, this is combination suitable for using \code{\link{multilevel.reliability}} }
\item{key}{A key vector showing which items have been reversed}
\item{ml}{The multilevel output} 
}
\references{
Cattell, R. B. (1964). Validity and reliability: A proposed more basic set of concepts. Journal of Educational Psychology, 55(1), 1 - 22. doi: 10.1037/h0046462

Cranford, J. A., Shrout, P. E., Iida, M., Rafaeli, E., Yip, T., \& Bolger, N. (2006). A procedure for evaluating sensitivity to within-person change: Can mood measures in diary studies detect change reliably? Personality and Social Psychology Bulletin, 32(7), 917-929.

DeSimone, J. A. (2015). New techniques for evaluating temporal consistency. Organizational Research Methods, 18(1), 133-152. doi: 10.1177/1094428114553061

Revelle, W. and Condon, D. Reliability  (in prep)

Revelle, W.  (in preparation) An introduction to psychometric theory with applications in {R}. Springer.  (Available online at \url{https://personality-project.org/r/book}). 

Shrout, P. E., & Lane, S. P. (2012). Psychometrics. In Handbook of research methods for studying daily life. Guilford Press.

Wood, D., Harms, P. D., Lowman, G. H., & DeSimone, J. A. (2017). Response speed and response consistency as mutually validating indicators of data quality in online samples. Social Psychological and Personality Science, 8(4), 454-464. doi: 10.1177/1948550617703168

}
\author{
William Revelle
}
\note{
lmer=TRUE is the default and will do the variance decomposition using lmer.  This will take some time.  For 3032 cases with 10 items from the msqR and sai data set, this takes 92 seconds, but just .63 seconds if lmer = FALSE.
For the 1895 subjects with repeated measures on the \code{\link[psychTools]{sai}}, it takes 85 seconds with lmer and .38 without out lmer.

In the case of just two tests (no items specified), the item based statistics (alpha, rqq, item.stats, scores, xy.df) are not reported. 

Two examples are given.  The first takes 200 cases from the \code{\link[psychTools]{sai}} data set.  Subjects were given the \code{link[psychTools]{sai}} twice with an intervening mood manipulation (four types of short film clips, with or without placebo/caffeine).  The test retest stability of the sai are based upon the 20 sai items.  The second example compares the scores of the 10 sai items that overlap with 10 items from the \code{\link[psychTools]{msqR}} data set from the same study.  \code{link[psychTools]{sai}} and  \code{\link[psychTools]{msqR}} were given immediately after each other and although the format differs slightly, can be seen as measures of dependability.  

}


\seealso{
 \code{\link{alpha}}, \code{\link{omega}} \code{\link{scoreItems}}, \code{\link{cor2}}
}
\examples{
\donttest{  #for faster compiling, dont test 
#lmer set to FALSE for speed.
#set lmer to TRUE to get variance components
sai.xray <- subset(psychTools::sai,psychTools::sai$study=="XRAY")
#The case where the two measures are identified by time
#automatically reverses items but throws a warning
stability <- testRetest(sai.xray[-c(1,3)],lmer=FALSE) 
stability  #show the results
#get a second data set
sai.xray1 <- subset(sai.xray,sai.xray$time==1)
msq.xray <- subset(psychTools::msqR,
 (psychTools::msqR$study=="XRAY") & (psychTools::msqR$time==1))
select <- colnames(sai.xray1)[is.element(colnames(sai.xray1 ),colnames(psychTools::msqR))] 

select <-select[-c(1:3)]  #get rid of the id information
#The case where the two times are in the form x, y

dependability <-  testRetest(sai.xray1,msq.xray,keys=select,lmer=FALSE)
dependability  #show the results
}
#now examine the Impulsivity subscale of the EPI
#use the epiR data set which includes epi.keys
data("epiR",package="psychTools")
#Imp <- selectFromKeys(epi.keys$Imp)   #fixed temporarily with 
Imp <- c("V1", "V3", "V8", "V10","V13" ,"V22", "V39" , "V5" , "V41")
imp.analysis <- testRetest(psychTools::epiR,select=Imp) #test-retest = .7, alpha=.51,.51 
imp.analysis
}
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
\keyword{ models }% __ONLY ONE__ keyword per line
