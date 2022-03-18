\name{anova.psych}
\alias{anova.psych}

\title{Model comparison for regression, mediation, and factor analysis}
\description{
When doing regressions from the data or from a correlation matrix using \code{\link{setCor}} or doing a mediation analysis using \code{link{mediate}}, it is useful to compare alternative models.  Since these are both regression models, the appropriate test is an Analysis of Variance.  Similar tests, using Chi Square may be done for factor analytic models. 
}
\usage{
\method{anova}{psych}(object,...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object from \code{\link{setCor}},  \code{\link{mediate}}, \code{\link{omega}}, or \code{\link{fa}}.
}
  \item{\dots}{More objects of the same type may be supplied here}

}
\details{
\code{\link{setCor}} returns the SE.residual and degrees of freedom.  These are converted to SSR and then an analysis of variance is used to compare two (or more) models. For \code{\link{omega}} or \code{\link{fa}} the change in the ML chisquare statistic as a function of change in df is reported.
}
\value{An ANOVA table comparing the models.}



\author{
Wiliam Revelle
}
\note{
The code has been adapted from the anova.lm function in stats and the anova.sem by John Fox.
}


\seealso{\code{\link{setCor}}, \code{\link{mediate}}, \code{\link{omega}}, \code{\link{fa}}

}
\examples{
if(require("psychTools")) {
m1 <- setCor(reaction ~ import, data = Tal_Or,std=FALSE)
m2 <- setCor(reaction ~ import+pmi, data = Tal_Or,std=FALSE)
m3 <- setCor(reaction ~ import+pmi + cond, data = Tal_Or,std=FALSE)
anova(m1,m2,m3)
}


#Several interesting test cases are taken from analyses of the Spengler data set
#Although the sample sizes are actually very large in the first wave,  I use the
#sample sizes from the last wave 
#This data set is actually in psychTools but is copied here until we can update psychTools
#We set the n.iter to be 50 instead of the default value of 5,000
if(require("psychTools")) {

 mod1 <- mediate(Income.50 ~ IQ + Parental+ (Ed.11) ,data=Spengler,
    n.obs = 1952, n.iter=50)
 mod2 <- mediate(Income.50 ~ IQ + Parental+ (Ed.11)  + (Income.11)
  ,data=Spengler,n.obs = 1952, n.iter=50)

#Now, compare these models
anova(mod1,mod2)
}

f3 <- fa(Thurstone,3,n.obs=213)  #we need to specifiy the n.obs for the test to work
f2 <- fa(Thurstone,2, n.obs=213)
anova(f2,f3)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line
