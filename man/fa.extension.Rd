\name{fa.extension}
\alias{fa.extension}
\alias{fa.extend}
\alias{faRegression}
\alias{faReg}
\title{Apply Dwyer's factor extension to find factor loadings for extended variables}

\description{Dwyer (1937) introduced a method for finding factor loadings for variables not included in the original analysis.  This is basically finding the unattenuated correlation of the extension variables with the factor scores.  An alternative, which does not correct for factor reliability was proposed by Gorsuch (1997). Both options are an application of exploratory factor analysis with extensions to new variables. Also useful for finding the validities of variables in the factor space.  
}
\usage{
fa.extension(Roe,fo,correct=TRUE)
fa.extend(r,nfactors=1,ov=NULL,ev=NULL,n.obs = NA, np.obs=NULL,
  correct=TRUE,rotate="oblimin",SMC=TRUE,   warnings=TRUE, fm="minres",
  alpha=.1,omega=FALSE,cor="cor",use="pairwise",cor.correct=.5,weight=NULL,
    missing=FALSE, smooth=TRUE, ...) 
faRegression(r,nfactors=1,ov=NULL,dv=NULL, n.obs = NA
 , np.obs=NULL,correct=TRUE,rotate="oblimin",SMC=TRUE,warnings=TRUE, fm="minres",alpha=.1,
 omega=FALSE,cor="cor",use="pairwise",cor.correct=.5,weight=NULL,smooth=TRUE, ...)
#this is just an alias for 
faReg(r,nfactors=1,ov=NULL,dv=NULL, n.obs = NA
 , np.obs=NULL,correct=TRUE,rotate="oblimin",SMC=TRUE,warnings=TRUE, fm="minres",alpha=.1,
 omega=FALSE,cor="cor",use="pairwise",cor.correct=.5,weight=NULL,smooth=TRUE, ...)
}

\arguments{
  \item{Roe}{The correlations of the original variables with the extended variables}
  \item{fo}{The output from the \code{\link{fa}} or \code{\link{omega}} functions applied to the original variables.}
  \item{correct}{correct=TRUE produces Dwyer's solution, correct=FALSE does not correct for factor reliability.  This is not quite the Gorsuch technique.}
\item{r}{A correlation or data matrix with all of the variables to be analyzed by fa.extend}
\item{ov}{The original variables to factor}
\item{ev}{The extension variables}
\item{dv}{The dependent variables if doing faRegression}
 \item{nfactors}{ Number of factors to extract, default is 1 }
  \item{n.obs}{Number of observations used to find the correlation matrix if using a correlation matrix.  Used for finding the goodness of fit statistics. Must be specified if using a correlaton matrix and finding confidence intervals.}
\item{np.obs}{Pairwise number of observations.  Required if using fm="minchi", suggested in other cases to estimate the empirical goodness of fit.}
   \item{rotate}{"none", "varimax", "quartimax",  "bentlerT",  "geominT" and "bifactor" are orthogonal rotations.  "promax", "oblimin", "simplimax", "bentlerQ,  "geominQ" and "biquartimin" and "cluster" are possible rotations or transformations of the solution.  The default is to do a oblimin transformation, although  versions prior to 2009 defaulted to varimax.}
  
  \item{SMC}{Use squared multiple correlations (SMC=TRUE) or use 1  as initial communality estimate. Try using 1 if imaginary eigen values are reported. If SMC is a vector of length the number of variables, then these values are used as starting values in the case of fm='pa'. }

  \item{warnings}{warnings=TRUE => warn if number of factors is too many }
  \item{fm}{factoring method  fm="minres" will do a minimum residual (OLS), fm="wls" will do a weighted least squares (WLS) solution, fm="gls" does a generalized weighted least squares (GLS), fm="pa" will do the principal factor solution, fm="ml" will do a maximum likelihood factor analysis. fm="minchi" will minimize the sample size weighted chi square when treating pairwise correlations with different number of subjects per pair.}
\item{alpha}{alpha level for the confidence intervals for RMSEA}
\item{omega}{Do the extension analysis for an omega type analysis}
\item{cor}{Pass the kind of correlation to fa (defaults to Pearson, can use mixed)}
\item{use}{Option for the cor function on how to handle missing data.}
\item{cor.correct}{The correction to be passed to mixed, tet, or polycor (defaults to .5)}
\item{weight}{Should we weight the variables? (see \code{\link{fa}})} 
\item{missing}{When finding factor scores, are missing data allowed?}
\item{smooth}{Should we smooth the correlation matrix -- smoothing produces bad output if there are NAs in the R matrix }
 \item{...}{Additional parameters, specifically, keys may be passed if using the target rotation, or delta if using geominQ, or whether to normalize if using Varimax}

}
\details{It is sometimes the case that factors are derived from a set of variables (the Fo factor loadings) and we want to see what the loadings of an extended set of variables (Fe) would be. Given the original correlation matrix Ro and the correlation of these original variables with the extension variables of Roe, it is a straight forward calculation to find the loadings Fe of the extended variables on the original factors.  This technique was developed by Dwyer (1937) for the case of adding new variables to a factor analysis without doing all the work over again. But, as discussed by Horn (1973) factor extension is also appropriate when one does not want to include the extension variables in the original factor analysis, but does want to see what the loadings would be anyway.

This could be done by estimating the factor scores and then finding the covariances of the extension variables with the factor scores. If raw data are available, this is done by \code{\link{faReg}}. But if the original data are not available, but just the covariance or correlation matrix is, then the use of \code{\link{fa.extension}} is most appropriate to estimate those correlations before doing the regression.  

The factor analysis results from either \code{\link{fa}} or \code{\link{omega}} functions applied to the original correlation matrix is extended to the extended variables given the correlations (Roe) of the extended variables with the original variables.

\code{\link{fa.extension}} assumes that the original factor solution was found by the \code{\link{fa}} function.

For a very nice discussion of the relationship between factor scores, correlation matrices, and the factor loadings in a factor extension, see Horn (1973).

The \code{\link{fa.extend}} function may be thought of as a "seeded" factor analysis.  That is, the variables in the original set are factored, this solution is then extended to the extension set, and the resulting output is presented as if both the original and extended variables were factored together.  This may also be done for an omega analysis. 

The example of  \code{\link{fa.extend}} compares the extended solution to a direct solution of all of the variables using \code{\link{factor.congruence}}. 

Another important use of factor extension is when parts of a correlation matrix are missing.  Thus suppose you have a correlation matrix formed of variables 1 .. i, j..n  where the first i variables correlate with each other (R11) and with the  j... n (R12) variables, but the elements of of the R22 (j..n) are missing.  That is, X has a missing data structure but we can find the correlation matrix R[1..n,1..n]  =cor(X,na.rm=TRUE) 

R[1..n,1..n] = \tabular{rr}{
   R[1..i,1..i] \tab R[1..i,j..n]  \cr
   R[j..n,1..i] \tab NA  \cr 
}
Factors can be found for the R11 matrix and then extended to the variables in the entire matrix. This allows for irt approaches to be applied even with significantly missing data.

f <- fa.extend(R,nf=1, ov=1:i,ev=j:n)

Then the loadings of n loadings of f may be found. 

Combining this procedure with \code{\link{fa2irt}} allows us to find the irt based parameters of the n variables, even though we have substantially incomplete data. 


Using the idea behind factor extension it is straightforward to apply these techniques to multiple regression, because the extended loadings are functionally beta weights.  . \code{\link{faRegression}} just organizes the extension analysis in terms  of a regression analysis.  If the raw data are available, it first finds the factor scores and then correlates these with the dependent variable in the standard regression approach.  But, if just the correlation matrix is available, it estimates the factor by dependent variable correlations by the use of \code{link{fa.extend}} before doing the regression.
}
\value{Factor Loadings of the exended variables on the original factors}

\references{Paul S. Dwyer (1937) The determination of the factor loadings of a given test from the known factor loadings of other tests. Psychometrika, 3, 173-178

 Gorsuch, Richard L. (1997) New procedure for extension analysis in exploratory factor analysis,  Educational and Psychological Measurement, 57, 725-740
 
 Horn, John L. (1973) On extension analysis and its relation to correlations between variables and factor scores.  Multivariate Behavioral Research, 8, (4), 477-489.
}
\author{William Revelle
}




\seealso{ See Also as \code{\link{fa}}, \code{\link{principal}}, \code{\link{Dwyer}}
}
\examples{ #The Dwyer Example
Ro <- Dwyer[1:7,1:7]
Roe <- Dwyer[1:7,8]
fo <- fa(Ro,2,rotate="none")
fe <- fa.extension(Roe,fo)

#an example from simulated data
set.seed(42) 
 d <- sim.item(12)    #two orthogonal factors 
 R <- cor(d)
 Ro <- R[c(1,2,4,5,7,8,10,11),c(1,2,4,5,7,8,10,11)]
 Roe <- R[c(1,2,4,5,7,8,10,11),c(3,6,9,12)]
 fo <- fa(Ro,2)
 fe <- fa.extension(Roe,fo)
 fa.diagram(fo,fe=fe)
 
 #alternatively just specify the original variables and the extension variables
 fe = fa.extend(R, 2, ov =c(1,2,4,5,7,8,10,11), ev=c(3,6,9,12))
 fa.diagram(fe$fo, fe = fe$fe)
 
 #create two correlated factors
 fx <- matrix(c(.9,.8,.7,.85,.75,.65,rep(0,12),.9,.8,.7,.85,.75,.65),ncol=2)
 Phi <- matrix(c(1,.6,.6,1),2)
 sim.data <- sim.structure(fx,Phi,n=1000,raw=TRUE)
 R <- cor(sim.data$observed)
 Ro <- R[c(1,2,4,5,7,8,10,11),c(1,2,4,5,7,8,10,11)]
 Roe <- R[c(1,2,4,5,7,8,10,11),c(3,6,9,12)]
 fo <- fa(Ro,2)
 fe <- fa.extension(Roe,fo)
 fa.diagram(fo,fe=fe)
 
 #now show how fa.extend works with the same data set
 #note that we have to make sure that the variables are in the order to do the factor congruence
 fe2 <- fa.extend(sim.data$observed,2,ov=c(1,2,4,5,7,8,10,11),ev=c(3,6,9,12))
 fa.diagram(fe2,main="factor analysis with extension variables")
 fa2 <- fa(sim.data$observed[,c(1,2,4,5,7,8,10,11,3,6,9,12)],2)
 factor.congruence(fe2,fa2)
 summary(fe2)
 
 #an example of extending an omega analysis
 
 
fload <- matrix(c(c(c(.9,.8,.7,.6),rep(0,20)),c(c(.9,.8,.7,.6),rep(0,20)),c(c(.9,.8,.7,.6),
        rep(0,20)),c(c(c(.9,.8,.7,.6),rep(0,20)),c(.9,.8,.7,.6))),ncol=5)
 gload <- matrix(rep(.7,5))
 five.factor <- sim.hierarchical(gload,fload,500,TRUE) #create sample data set
 ss <- c(1,2,3,5,6,7,9,10,11,13,14,15,17,18,19)
 Ro <- cor(five.factor$observed[,ss])
 Re <- cor(five.factor$observed[,ss],five.factor$observed[,-ss])
 om5 <-omega(Ro,5)   #the omega analysis
 om.extend <- fa.extension(Re,om5) #the extension analysis
 om.extend #show it
 #now, include it in an omega diagram
 combined.om <- rbind(om5$schmid$sl[,1:ncol(om.extend$loadings)],om.extend$loadings)
 class(combined.om) <-c("psych","extend")
 omega.diagram(combined.om,main="Extended Omega") 
 
 #show how to use fa.extend to do regression analyses with the raw data
 b5 <- faReg (bfi, nfactors = 5, ov =1:25, dv =26:28)
 extension.diagram(b5)
 R <-cor(bfi,use="pairwise")
 b5.r <- faReg(R, nfactors = 5, ov =1:25, dv =26:28) # not identical to b5
 round(b5$regression$coefficients - b5.r$regression$coefficients,2)
}
\keyword{ multivariate }