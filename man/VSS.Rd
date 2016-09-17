\name{VSS}
\alias{vss}
\alias{VSS}
\alias{MAP}
\alias{nfactors}

\title{ Apply the Very Simple Structure, MAP, and other criteria to determine the appropriate number of factors.}
\description{There are multiple ways to determine the appropriate number of factors in exploratory factor analysis. Routines for the Very Simple Structure (VSS) criterion allow one to compare solutions of varying complexity and for different number of factors. Graphic output indicates the "optimal" number of factors for different levels of complexity.  The Velicer MAP criterion is another good choice. \code{\link{nfactors}} finds and plots several of these alternative estimates.
}
\usage{
vss(x, n = 8, rotate = "varimax", diagonal = FALSE, fm = "minres", 
n.obs=NULL,plot=TRUE,title="Very Simple Structure",use="pairwise",cor="cor",...)
VSS(x, n = 8, rotate = "varimax", diagonal = FALSE, fm = "minres", 
n.obs=NULL,plot=TRUE,title="Very Simple Structure",use="pairwise",cor="cor",...)
nfactors(x,n=20,rotate="varimax",diagonal=FALSE,fm="minres",n.obs=NULL,
             title="Number of Factors",pch=16,use="pairwise", cor="cor",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a correlation matrix or a data matrix}
  \item{n}{Number of factors to extract -- should be more than hypothesized! }
  \item{rotate}{ what rotation to use c("none", "varimax",  "oblimin","promax")}
  \item{diagonal}{Should we fit the diagonal as well }
  \item{fm}{factoring method -- fm="pa"  Principal Axis Factor Analysis, fm = "minres" minimum residual (OLS) factoring fm="mle"  Maximum Likelihood FA, fm="pc" Principal Components" }
  \item{n.obs}{Number of observations if doing a factor analysis of correlation matrix.  This value is ignored by VSS but is necessary for the ML factor analysis package.}
  \item{plot}{plot=TRUE  Automatically call VSS.plot with the VSS output, otherwise don't plot}
\item{title}{a title to be passed on to VSS.plot}
\item{pch}{the plot character for the nfactors plots}
\item{use}{If doing covariances or Pearson R, should we use "pairwise" or "complete cases"}
\item{cor}{What kind of correlation to find, defaults to Pearson but see fa for the choices}
\item{\dots}{parameters to pass to the factor analysis program 
      The most important of these is if using a correlation matrix is covmat= xx}
}
\details{Determining the most interpretable number of factors from a factor analysis is perhaps one of the greatest challenges in factor analysis.  There are many solutions to this problem, none of which is uniformly the best.  "Solving the number of factors problem is easy, I do it everyday before breakfast."  But knowing the right solution is harder. (Horn and Engstrom, 1979) (Henry Kaiser in personal communication with J.L. Horn, as cited by Horn and Engstrom, 1979, MBR p 283).  

Techniques most commonly used include

	
1)  Extracting factors until the chi square of the residual matrix is not significant.

2) Extracting factors until the change in chi square from factor n to factor n+1 is not significant.

3) Extracting factors until the eigen values of the real data are less than the corresponding eigen values of a random data set of the same size (parallel analysis) \code{\link{fa.parallel}}.

4) Plotting the magnitude of the successive eigen values and applying the scree test (a sudden drop in eigen values analogous to the change in slope seen when scrambling up the talus slope of a mountain and approaching the rock face.

5) Extracting principal components until the eigen value < 1. 

6) Extracting factors as long as they are interpetable.
	
7) Using the Very Simple Structure Criterion (VSS).

8) Using Wayne Velicer's Minimum Average Partial (MAP) criterion. 
	

 Each of the procedures has its advantages and disadvantages.  Using either the chi square test or the change in square test is, of course, sensitive to the number of subjects and leads to the nonsensical condition that if one wants to find many factors, one simply runs more subjects. Parallel analysis is partially sensitive to sample size in that for large samples the eigen values of random factors will be very small.  The scree test is quite appealling but can lead to differences of interpretation as to when the scree "breaks". The eigen value of 1 rule, although the default for many programs, seems to be a rough way of dividing the number of variables by 3.  Extracting interpretable factors means that the number of factors reflects the investigators creativity more than the data.  VSS, while very simple to understand, will not work very well if the data are very factorially complex. (Simulations suggests it will work fine if the complexities of some of the items are no more than 2).



Most users of factor analysis tend to interpret factor output by focusing their attention on the largest loadings for every variable and ignoring the smaller ones.  Very Simple Structure operationalizes this tendency by  comparing the original correlation matrix to that reproduced by a simplified version (S) of the original factor matrix (F).  R = SS' + U2.   S is composed of just the c greatest (in absolute value) loadings for each variable.  C (or complexity) is a parameter of the model and may vary from 1 to the number of factors.  

The VSS criterion compares the fit of the simplified model to the original correlations: VSS = 1 -sumsquares(r*)/sumsquares(r)  where R* is the residual matrix R* = R - SS' and r* and r are the elements of R* and R respectively. 

VSS for a given complexity will tend to peak at the optimal (most interpretable) number of factors (Revelle and Rocklin, 1979). 

Although originally written in Fortran for main frame computers, VSS has been adapted to micro computers (e.g., Macintosh OS 6-9) using Pascal. We now release R code for calculating VSS. 

Note that if using a correlation matrix (e.g., my.matrix) and doing a factor analysis, the parameters n.obs should be specified for the factor analysis:
e.g., the call is VSS(my.matrix,n.obs=500).  Otherwise it defaults to 1000. 

Wayne Velicer's MAP criterion has been added as an additional test for the optimal number of components to extract.  Note that VSS and MAP will not always agree as to the optimal number.

The nfactors function will do a VSS, find MAP, and report a number of other criteria (e.g., BIC, complexity, chi square, ...)



A variety of rotation options are available. These include varimax, promax, and oblimin. Others can be added.  Suggestions are welcome.
}
\value{ A data.frame with entries:
 map:  Velicer's MAP values (lower values are better) \cr
 dof:  degrees of freedom    (if using FA) \cr
 chisq: chi square (from the factor analysis output  (if using FA) \cr
 prob: probability of residual matrix > 0  (if using FA) \cr
 sqresid: squared residual correlations\cr
 RMSEA: the RMSEA for each number of factors \cr
 BIC: the BIC for each number of factors \cr
 eChiSq: the empirically found chi square \cr
 eRMS: Empirically found mean residual \cr
 eCRMS: Empirically found mean residual corrected for df \cr
 eBIC: The empirically found BIC based upon the eChiSq \cr
 fit: factor fit of the complete model\cr
 cfit.1: VSS fit of complexity 1\cr
 cfit.2: VSS fit of complexity 2 \cr
 ... \cr
 cfit.8: VSS fit of complexity 8\cr
 cresidiual.1: sum squared residual correlations for complexity 1\cr
 ...:  sum squared residual correlations for complexity 2 ..8\cr
 
}
\references{ \url{http://personality-project.org/r/vss.html},
Revelle, W. An introduction to psychometric theory with applications in R (in prep) Springer. Draft chapters available at  \url{http://personality-project.org/r/book/}

Revelle, W. and Rocklin, T. 1979, Very Simple Structure: an Alternative Procedure for Estimating the Optimal Number of Interpretable Factors, Multivariate Behavioral Research, 14, 403-414. 
\url{http://personality-project.org/revelle/publications/vss.pdf}

Velicer, W. (1976) Determining the number of components from the matrix of partial correlations. Psychometrika, 41, 321-327.
	
}

\author{William Revelle}

\seealso{  \code{\link{VSS.plot}},  \code{\link{ICLUST}}, \code{\link{omega}},
 \code{\link{fa.parallel}}}
\examples{

#test.data <- Harman74.cor$cov
#my.vss <- VSS(test.data,title="VSS of 24 mental tests")      
#print(my.vss[,1:12],digits =2) 
#VSS.plot(my.vss, title="VSS of 24 mental tests")

#now, some simulated data with two factors
#VSS(sim.circ(nvar=24),fm="minres" ,title="VSS of 24 circumplex variables")
VSS(sim.item(nvar=24),fm="minres" ,title="VSS of 24 simple structure variables")
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
