\name{esem}
\alias{esem}
\alias{esem.diagram}
\alias{interbattery}

\title{Perform and Exploratory Structural Equation Model (ESEM) by using factor extension techniques}

\description{Structural Equation Modeling (SEM) is a powerful tool for confirming multivariate structures and is well done by the lavaan, sem, or OpenMx packages. Because they are confirmatory, SEM  models test specific models.  Exploratory Structural Equation Modeling (ESEM), on the other hand, takes a more exploratory approach.  By using factor extension, it is possible to extend the factors of one set of variables (X) into the variable space of another set (Y). Using this technique, it is then possible to estimate the correlations between the two sets of latent variables, much the way normal SEM would do.  Based upon exploratory factor analysis (EFA) this approach provides a quick and easy approach to do exploratory structural equation modeling.  
}
\usage{
esem(r, varsX, varsY, nfX = 1, nfY = 1, n.obs = NULL, fm = "minres", 
      rotate = "oblimin", plot = TRUE, cor = "cor", use = "pairwise",weight=NULL, ...)
esem.diagram(esem=NULL,labels=NULL,cut=.3,errors=FALSE,simple=TRUE,
regression=FALSE,lr=TRUE, digits=1,e.size=.1,adj=2,
     main="Exploratory Structural Model", ...)
interbattery(r, varsX, varsY, nfX = 1, nfY = 1, n.obs = NULL,cor = "cor", 
       use = "pairwise",weight=NULL)
}
\arguments{
 \item{r}{A correlation matrix or a raw data matrix suitable for factor analysis}
\item{varsX}{The variables defining set X}
\item{varsY}{The variables defining set Y}
 \item{nfX}{The number of factors to extract for the X variables}
 \item{nfY}{The number of factors to extract for the Y variables}
  \item{n.obs}{Number of observations (needed for eBIC and chi square), can be ignored.}
  \item{fm}{The factor  method to use, e.g., "minres", "mle" etc.  (see fa for details)}
  \item{rotate}{Which rotation to use. (see fa for details)}
  \item{plot}{If TRUE, draw the esem.diagram}
  \item{cor}{What options for to use for correlations (see fa for details)}
  \item{use}{"pairwise" for pairwise complete data, for other options see cor}
  \item{weight}{Weights to apply to cases when finding wt.cov}
 \item{\dots}{other parameters to pass to fa or to esem.diagram functions.}
 \item{esem}{The object returned from esem and passed to esem.diagram}
   \item{labels}{ Variable labels }
  \item{cut}{ Loadings with abs(loading) > cut will be shown }
  \item{simple}{Only the biggest loading per item is shown}
   \item{errors}{include error estimates (as arrows)}
  \item{e.size}{size of ellipses (adjusted by the number of variables)}
  \item{digits}{Round coefficient to digits}
  \item{adj}{loadings are adjusted by factor number mod adj to decrease likelihood of overlap}
   \item{main}{ Graphic title, defaults to "Exploratory Structural Model" }
   \item{lr}{draw the graphic left to right (TRUE) or top to bottom (FALSE)}
   \item{regression}{Not yet implemented}
   }
\details{
Factor analysis as implemented in \code{\link{fa}} attempts to  summarize the covariance (correlational) structure of a set of variables with a small set of latent variables or ``factors".  This solution may be `extended' into a larger space with more variables without changing the original solution (see \code{\link{fa.extension}}.  Similarly,  the factors of a second set of variables  (the Y set) may be extended into the original (X ) set.  Doing so allows two independent measurement models, a measurement model for X and a measurement model for Y.  These two sets of latent variables may then be correlated  for an Exploratory Structural Equation Model.  (This is exploratory because it is based upon exploratory factor analysis (EFA) rather than a confirmatory factor model (CFA) using more traditional Structural Equation Modeling packages such as sem, lavaan, or Mx.)

Although the output seems very similar to that of a normal EFA using  \code{\link{fa}}, it is actually two independent factor analyses (of the X and the Y sets) that are then mutually extended into each other.  That is, the loadings and structure matrices from sets X and Y are merely combined, and the correlations between the two sets of factors are found.

Interbattery factor analysis was developed by Tucker (1958) as a way of comparing the factors in common to two batteries of tests.   (Currently under development and not yet complete). Using some straight forward linear algebra It is easy to find the factors of the intercorrelations between the two sets of variables.  This does not require estimating communalities and is highly related to the procedures of canonical correlation.  

The difference between the esem and the interbattery approach is that the first factors the X set and then relates those factors to factors of the Y set.  Interbattery factor analysis, on the other hand, tries to find one set of factors that links both sets but is still distinct from factoring both sets together.  
}
\value{
\item{communality}{The amount of variance in each of the X and Y variables accounted for by the total model.}
\item{sumsq}{The amount of variance accounted for by each factor -- independent of the other factors.}
\item{dof}{Degrees of freedom of the model}
\item{null.dof}{Degrees of freedom of the null  model (the correlation matrix)}
\item{ENull}{chi square of the null model}
\item{chi}{chi square of the model.  This is found by examining the size of the residuals compared to their standard error.}
\item{rms}{The root mean square of the residuals.}
\item{nh}{Harmonic sample size if using min.chi for factor extraction.}
\item{EPVAL}{Probability  of the Emprical Chi Square given the hypothesis of an identity matrix.}
\item{crms}{Adjusted root mean square residual}
\item{EBIC}{When normal theory fails (e.g., in the case of non-positive definite matrices), it useful to examine the empirically derived EBIC based upon the empirical \eqn{\chi^2}{chi^2} - 2 df. }
\item{ESABIC}{Sample size adjusted empirical BIC}
\item{fit}{sum of squared residuals versus sum of squared original values}
\item{fit.off}{fit applied to the off diagonal elements}
\item{sd}{standard deviation of the residuals}
\item{factors}{Number of factors extracted}
\item{complexity}{Item complexity}
\item{n.obs}{Number of total observations}
\item{loadings}{The factor pattern matrix for the combined X and Y factors}
\item{Structure}{The factor structure matrix for the combined X and Y factors} 
\item{loadsX}{Just the X set of loadings (pattern) without  the extension variables.}
\item{loadsY}{Just the Y set of loadings (pattern) without the extension variables.}      \item{PhiX}{The correlations of the X factors}
\item{PhiY}{the correlations of the Y factors}
\item{Phi}{the correlations of the X and Y factors within the selves and across sets.}
\item{fm}{The factor method used}
\item{fx}{The complete factor analysis output for the X set}
\item{fy}{The complete factor analysis output for the Y set}
\item{residual}{The residual correlation matrix (R - model).  May be examined by a call to residual().}    
\item{Call}{Echo back the original call to the function.}
\item{model}{model code for SEM and for lavaan to do subsequent confirmatory modeling}
}
\references{
Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 

Tucker, Ledyard (1958) An inter-battery method of factor analysis, Psychometrika, 23, 111-136.


}
\author{William Revelle}

\note{Developed September, 2016, revised December, 2018 to produce code for lavaan and sem.  the \code{\link{esem}} and \code{\link{esem.diagram}} functions. Suggestions or comments are most welcome. 
}


\seealso{ \code{\link{principal}} for principal components analysis (PCA).  PCA will give very similar solutions to factor analysis when there are many variables.  The differences become more salient as the number variables decrease.  The PCA and FA models are actually very different and should not be confused.  One is a model of the observed variables, the other is a model of latent variables.  

\code{\link{irt.fa}} for Item Response Theory analyses using factor analysis, using the two parameter IRT equivalent of loadings and difficulties.

\code{\link{VSS}} will produce the Very Simple Structure (VSS) and MAP criteria for the number of factors, \code{\link{nfactors}} to compare many different factor criteria.

\code{\link{ICLUST}} will do a hierarchical cluster analysis alternative to factor analysis or principal components analysis.

 \code{\link{predict.psych}} to find predicted scores based upon new data, \code{\link{fa.extension}} to extend the factor solution to new variables, \code{\link{omega}} for hierarchical factor analysis with one general factor. 
 code{\link{fa.multi}} for hierarchical factor analysis with an arbitrary number of higher order factors. 
 

\code{\link{fa.sort}} will sort the factor loadings into echelon form. \code{\link{fa.organize}} will reorganize the factor pattern matrix into any arbitrary order of factors and items.  

\code{\link{KMO}} and \code{\link{cortest.bartlett}} for various tests that some people like. 

\code{\link{factor2cluster}} will prepare unit weighted scoring keys of the factors that can be used with \code{\link{scoreItems}}.

\code{\link{fa.lookup}} will print the factor analysis loadings matrix along with the item ``content" taken from a dictionary of items.  This is useful when examining the meaning of the factors.  

\code{\link{anova.psych}} allows for testing the difference between two (presumably nested) factor models .

}
\examples{
#make up a sem like problem using sim.structure
fx <-matrix(c( .9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)  
fy <- matrix(c(.6,.5,.4),ncol=1)
rownames(fx) <- c("V","Q","A","nach","Anx")
rownames(fy)<- c("gpa","Pre","MA")
Phi <-matrix( c(1,0,.7,.0,1,.7,.7,.7,1),ncol=3)
gre.gpa <- sim.structural(fx,Phi,fy)
print(gre.gpa)

#now esem it:
example <- esem(gre.gpa$model,varsX=1:5,varsY=6:8,nfX=2,nfY=1,n.obs=1000,plot=FALSE)
example
esem.diagram(example,simple=FALSE)

#compare two alternative solutions to the first 2 factors of the neo.
#solution 1 is the normal 2 factor solution.
#solution 2 is an esem with 1 factor for the first 6 variables, and 1 for the second 6.
f2 <- fa(neo[1:12,1:12],2)
es2 <- esem(neo,1:6,7:12,1,1)
summary(f2)
summary(es2)
fa.congruence(f2,es2)

interbattery(Thurstone.9,1:4,5:9,2,2)  #compare to the solution of Tucker.  We are not there yet.

}

\keyword{ multivariate }
\keyword{ models}