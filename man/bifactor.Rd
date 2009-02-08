\name{bifactor}
\alias{bifactor}
\alias{Holzinger}
\alias{Reise}
\alias{Thurstone}
\docType{data}
\title{Three data sets showing a bifactor solution.}
\description{Holzinger-Swineford (1937) introduced the bifactor model of a general factor and uncorrelated group factors. The Holzinger correlation matrix is a 14 * 14 matrix from their paper. The Reise data set is 16 * 16 correlation matrix of mental health items. The Thurstone correlation matrix is a 9 * 9. 
}
\usage{data(bifactor)}

\details{Holzinger and Swineford introduced the bifactor model (one general factor and several group factors) for mental abilities.  This is a nice demonstration data set of a hierarchical factor structure that can be analyzed using the \code{\link{omega}} function or using sem. The bifactor model is typically used in measures of cognitive ability.

The 14 variables are ordered to reflect 3 spatial tests, 3 mental speed tests, 4 motor speed tests, and 4 verbal tests.  The sample size is 355.

More recent applications of the bifactor model are to the measurement of psychological status. The Reise data set is a correlation matrix based upon >35,000 observations to the Consumer Assessment of Health Care Provideers and Systems survey instrument. Reise, Morizot, and Hays (2007) describe a bifactor solution based upon 1,000 cases.    

The five factors from Reise et al. reflect Getting care quickly (1-3), Doctor communicates well (4-7), Courteous and helpful staff (8,9), Getting needed care (10-13), and and Health plan customer service (14-16). 

The 9 variable Thurstone problem is discussed in detail by R. P. McDonald (1999) and is used as example in the sem package.

\itemize{
\item Holzinger: 14 x 14 correlation matrix of ability tests, N = 355

\item Reise: 16 x 16 correlation matrix of health satisfaction items.  N = 35,000

\item Thurstone: 9 x 9 correlation matrix of ability tests, N = 213

}
}

\source{Holzinger:  Holzinger and Swineford (1937) \cr
Reise: Steve Reise (personal communication) \cr
sem help page (for Thurstone) 
}
\references{ Holzinger, Karl and Swineford, Frances (1937) The Bi-factor method.  Psychometrika, 2, 41-54

McDonald, Roderick P. (1999) Test theory: A unified treatment. L. Erlbaum Associates
Mahwah, N.J.

Reise, Steven and Morizot, Julien and Hays, Ron (2007) The role of the bifactor model in resolving dimensionality issues in health outcomes measures. Quality of Life Research. 16, 19-31.

}
\examples{
data(bifactor)
holz <- omega(Holzinger,4, title = "14 ability tests from Holzinger-Swineford")
bf <- omega(Reise,5,title="16 health items from Reise") 
omega(Reise,5,labels=colnames(Reise),title="16 health items from Reise")
thur.bf <- omega(Thurstone,title="9 variables from Thurstone")

}
\keyword{datasets}
