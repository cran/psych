\name{VSS}
\alias{VSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Apply the Very Simple Structure criterion to determine the appropriate number of factors.}
\description{There are multiple ways to determine the appropriate number of factors in exploratory factor analysis. Routines for the Very Simple Structure (VSS) criterion allow one to compare solutions of varying complexity and for different number of factors. Graphic output indicates the "optimal" number of factors for different levels of complexity. 
}
\usage{
VSS(x, n = 8, rotate = "varimax", diagonal = FALSE, pc = "pa", n.obs=1000,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a correlation matrix or a data matrix}
  \item{n}{Number of factors to extract -- should be more than hypothesized! }
  \item{rotate}{ what rotation to use c("none", "varimax", "promax")}
  \item{diagonal}{Should we fit the diagonal as well }
  \item{pc}{pc="pa"  Principal Axis Factor Analysis, pc="mle"  Maximum Likelihood FA, pc="pc" Principal Components" }
  \item{n.obs}{Number of observations if doing a factor analysis of correlation matrix.  This value is ignored by VSS but is necessary for the ML factor analysis package.}
  \item{\dots}{parameters to pass to the factor analysis program 
      The most important of these is if using a correlation matrix is covmat= xx}
}
\details{Determining the most interpretable number of factors from a factor analysis is perhaps one of the greatest challenges in factor analysis.  There are many solutions to this problem, none of which is uniformly the best.  "Solving the number of factors problem is easy, I do it everyday before breakfast.  But knowing the right solution is harder" (Kaiser, 195x).  

Techniques most commonly used include

	
1)  Extracting factors until the chi square of the residual matrix is not significant.

2) Extracting factors until the change in chi square from factor n to factor n+1 is not significant.

3) Extracting factors until the eigen values of the real data are less than the corresponding eigen values of a random data set of the same size (parallel analysis).

4) Plotting the magnitude of the successive eigen values and applying the scree test (a sudden drop in eigen values analogous to the change in slope seen when scrambling up the talus slope of a mountain and approaching the rock face.

5) Extracting principal components until the eigen value <1. <

6) Extracting factors as long as they are interpetable.
	
7) Using the Very Structure Criterion.
	

 Each of the procedures has its advantages and disadvantages.  Using either the chi square test or the change in square test is, of course, sensitive to the number of subjects and leads to the nonsensical condition that if one wants to find many factors, one simlpy runs more subjects. Parallel analysis is partially sensitive to sample size in that for large samples the eigen values of random factors will be very small.  The scree test is quite appealling but can lead to differences of interpretation as to when the scree "breaks". The eigen value of 1 rule, although the default for many programs, seems to be a rough way of dividing the number of variables by 3.  Extracting interpretable factors means that the number of factors reflects the investigators creativity more than the data.  VSS, while very simple to understand, will not work very well if the data are very factorially complex. (Simulations suggests it will work fine if the complexities of some of the items are no more than 2).



Most users of factor analysis tend to interpret factor output by focusing their attention on the largest loadings for every variable and ignoring the smaller ones.  Very Simple Structure operationalizes this tendency by  comparing the original correlation matrix to that reproduced by a simplified version (S) of the original factor matrix (F).  R = SS' + U2.   S is composed of just the c greatest (in absolute value) loadings for each variable.  C (or complexity) is a parameter of the model and may vary from 1 to the number of factors.  

The VSS criterion compares the fit of the simplified model to the original correlations: VSS = 1 -sumsquares(r*)/sumsquares(r)  where R* is the residual matrix R* = R - SS' and r* and r are the elements of R* and R respectively. 

VSS for a given complexity will tend to peak at the optimal (most interpretable) number of factors (Revelle and Rocklin, 1979). 

Although originally written in Fortran for main frame computers, VSS has been adapted to micro computers (e.g., Macintosh OS 6-9) using Pascal. We now release R code for calculating VSS. 

Note that if using a correlation matrix (e.g., my.matrix) and doing a factor analysis, the parameters n.obs should be specified for the factor analysis:
the call is VSS(my.matrix,n.obs=500).  Otherwise it defaults to 1000. 


}
\value{ A data.frame with entries:
 dof:  degrees of freedom    (if using FA) \cr
 chisq: chi square (from the factor analysis output  (if using FA) \cr
 prob: probability of residual matrix > 0  (if using FA) \cr
 sqresid: squared residual correlations\cr
 fit: factor fit of the complete model\cr
 cfit.1: VSS fit of complexity 1\cr
 cfit.2: VSS fit of complexity 2 \cr
 ... \cr
 cfit.8: VSS fit of complexity 8\cr
 cresidiual.1: sum squared residual correlations for complexity 1\cr
 ...:  sum squared residual correlations for complexity 2 ..8\cr
 
}
\references{ \url{http://personality-project.org/r/vss.html}
see also Revelle, W. and Rocklin, T. 1979, Very Simple Structure: an Alternative Procedure for Estimating the Optimal Number of Interpretable Factors, Multivariate Behavioral Research, 14, 403-414. \url{http://personality-project.org/revelle/publications/vss.pdf}}

\author{William Revelle \cr
        Department of Psychology \cr
        Northwestern University \cr
        Evanston, Illiniois \cr


Maintainer: William Revelle <revelle@northwestern.edu>

}

\seealso{  \code{\link{VSS.plot}},  \code{\link{ICLUST}}, \code{\link{omega}}}
\examples{
\dontrun{
test.data <- Harman74.cor$cov
my.vss <- VSS(test.data)         #suggests that 4 factor complexity two solution is optimal
print(my.vss[,1:12],digits =2) 
VSS.plot(my.vss)                 #see graphic window for a plot

#produces this output
#  dof chisq     prob sqresid  fit cfit.1 cfit.2 cfit.3 cfit.4 cfit.5 cfit.6 cfit.7
#1 252  4583  0.0e+00    17.2 0.79   0.79   0.00   0.00   0.00   0.00   0.00   0.00
#2 229  3105  0.0e+00    12.9 0.84   0.75   0.84   0.00   0.00   0.00   0.00   0.00
#4 186  1689 2.3e-240     8.0 0.90   0.66   0.87   0.90   0.90   0.00   0.00   0.00
#5 166  1398 9.3e-194     7.3 0.91   0.68   0.86   0.90   0.91   0.91   0.00   0.00
#6 147  1183 2.9e-161     6.5 0.92   0.53   0.83   0.88   0.91   0.92   0.92   0.00
#7 129  1002 5.8e-135     5.7 0.93   0.47   0.78   0.88   0.91   0.92   0.93   0.93
#8 112   803 5.3e-105     5.3 0.94   0.49   0.76   0.86   0.90   0.92   0.93   0.93

#compare the above solution to a "varimax" rotated solution which suggests 1 factor (g)

my.vss <- VSS(test.data,rotate="varimax")         #suggests that 1 factor complexity one solution is optimal
print(my.vss[,1:14],digits =2) 
VSS.plot(my.vss)                 #see graphic window for a plot

#  dof chisq     prob sqresid  fit cfit.1 cfit.2 cfit.3 cfit.4 cfit.5 cfit.6 cfit.7 cfit.8 cresidual.1
#1 252  4583  0.0e+00    17.2 0.79   0.79   0.00   0.00    0.0   0.00   0.00   0.00   0.00          17
#2 229  3105  0.0e+00    12.9 0.84   0.54   0.84   0.00    0.0   0.00   0.00   0.00   0.00          38
#3 207  2193  0.0e+00    10.1 0.88   0.46   0.79   0.88    0.0   0.00   0.00   0.00   0.00          45
#4 186  1689 2.3e-240     8.0 0.90   0.42   0.73   0.87    0.9   0.00   0.00   0.00   0.00          48
#5 166  1398 9.3e-194     7.3 0.91   0.40   0.70   0.86    0.9   0.91   0.00   0.00   0.00          50
#6 147  1183 2.9e-161     6.5 0.92   0.39   0.69   0.86    0.9   0.92   0.92   0.00   0.00          51
#7 129  1002 5.8e-135     5.7 0.93   0.39   0.70   0.84    0.9   0.92   0.93   0.93   0.00          50
#8 112   803 5.3e-105     5.3 0.94   0.39   0.69   0.83    0.9   0.92   0.93   0.93   0.94          50
}

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
