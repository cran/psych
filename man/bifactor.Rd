\name{bifactor}
\alias{bifactor}
\alias{Bechtoldt.1}
\alias{Bechtoldt.2}
\alias{Bechtoldt}
\alias{Holzinger}
\alias{Reise}
\alias{Thurstone}
\alias{Thurstone.33}
\docType{data}
\title{Six data sets showing a bifactor solution.}
\description{Holzinger-Swineford (1937) introduced the bifactor model of a general factor and uncorrelated group factors. The Holzinger correlation matrix is a 14 * 14 matrix from their paper. The Thurstone correlation matrix is a 9 * 9 matrix of correlations of ability items.  The Reise data set is 16 * 16 correlation matrix of mental health items. The Bechtholdt data sets are both 17 x 17 correlation matrices of ability tests.  
}
\usage{data(bifactor)}

\details{Holzinger and Swineford (1937) introduced the bifactor model (one general factor and several group factors) for mental abilities.  This is a nice demonstration data set of a hierarchical factor structure that can be analyzed using the \code{\link{omega}} function or using sem. The bifactor model is typically used in measures of cognitive ability.

The 14 variables are ordered to reflect 3 spatial tests, 3 mental speed tests, 4 motor speed tests, and 4 verbal tests.  The sample size is 355.

Another classic data set is the 9 variable Thurstone problem which is discussed in detail by R. P. McDonald (1985, 1999) and and is used as example in the sem package as well as in the PROC CALIS manual for SAS.  These nine tests were grouped by Thurstone, 1941 (based on other data) into three factors: Verbal Comprehension, Word Fluency, and Reasoning. The original data came from Thurstone and Thurstone (1941) but were reanalyzed by Bechthold (1961) who broke the data set into two. McDonald, in turn, selected these nine variables from a larger set of 17. The sample size is 213.

Another set of 9 cognitive variables attributed to Thurstone (1933) is the data set of 4,175 students reported by Professor Brigham of Princeton to the College Entrance Examination Board.  This set does not show a clear bifactor solution but is included as a demonstration of the differences between a maximimum likelihood factor analysis solution versus a principal axis factor solution.

More recent applications of the bifactor model are to the measurement of psychological status. The Reise data set is a correlation matrix based upon >35,000 observations to the Consumer Assessment of Health Care Provideers and Systems survey instrument. Reise, Morizot, and Hays (2007) describe a bifactor solution based upon 1,000 cases.    

The five factors from Reise et al. reflect Getting care quickly (1-3), Doctor communicates well (4-7), Courteous and helpful staff (8,9), Getting needed care (10-13), and and Health plan customer service (14-16). 

The two Bechtoldt data sets are two samples from Thurstone (1941).  They include 17 variables, 9 of which were used by McDonald to form the Thurstone data set. The sample sizes are 212 and 213 respectively. The six proposed factors reflect  memory, verbal, words,  space, number and reasoning with three markers for all expect the rote memory factor. 9 variables from this set appear in Thurstone.


\itemize{

\item Bechtoldt.1: 17 x 17 correlation matrix of ability tests, N = 212.

\item Bechtoldt.2: 17 x 17 correlation matrix of ability tests, N = 213.

\item Holzinger: 14 x 14 correlation matrix of ability tests, N = 355

\item Reise: 16 x 16 correlation matrix of health satisfaction items.  N = 35,000

\item Thurstone: 9 x 9 correlation matrix of ability tests, N = 213

\item Thurstone.33: Another 9 x 9 correlation matrix of ability items, N=4175

}
}

\source{Holzinger:  Holzinger and Swineford (1937) \cr
Reise: Steve Reise (personal communication) \cr
sem help page (for Thurstone) 
}
\references{ 
Bechtoldt, Harold, (1961). An empirical study of the factor analysis stability hypothesis. Psychometrika, 26, 405-432.

Holzinger, Karl and Swineford, Frances (1937) The Bi-factor method.  Psychometrika, 2, 41-54

McDonald, Roderick P. (1999) Test theory: A unified treatment. L. Erlbaum Associates. Mahwah, N.J.

Reise, Steven and Morizot, Julien and Hays, Ron (2007) The role of the bifactor model in resolving dimensionality issues in health outcomes measures. Quality of Life Research. 16, 19-31.

Thurstone, Louis Leon (1933) The theory of multiple factors. Edwards Brothers, Inc.  Ann Arbor

Thurstone, Louis Leon and Thurstone, Thelma (Gwinn). (1941) Factorial studies of intelligence. The University of Chicago Press. Chicago, Il.


}
\examples{
data(bifactor)
holz <- omega(Holzinger,4, title = "14 ability tests from Holzinger-Swineford")
bf <- omega(Reise,5,title="16 health items from Reise") 
omega(Reise,5,labels=colnames(Reise),title="16 health items from Reise")
thur.bf <- omega(Thurstone,title="9 variables from Thurstone")

}
\keyword{datasets}
