\name{tetrachoric}
\alias{tetrachoric}
\alias{tetrachor}
\alias{polychoric}
\alias{biserial}
\alias{polydi}
\alias{polyserial}
\alias{poly.mat}

\title{Tetrachoric, polychoric, biserial and polyserial correlations from various types of input}

\description{The tetrachoric correlation is the inferred Pearson Correlation from a two x two table with the assumption of bivariate normality. The polychoric correlation generalizes this to the n x m table. Particularly important when doing Item Response Theory or converting comorbidity statistics using normal theory to correlations. Input may be a 2 x 2 table of cell frequencies, a vector of cell frequencies, or a data.frame or matrix of dichotomous data (for tetrachoric) or of numeric data (for polychoric).
The biserial correlation is between a continuous y variable and a dichotmous x variable, which is assumed to have resulted from a dichotomized normal variable. Biserial is a special case of the polyserial correlation, which is the inferred latent correlation between a continuous variable (X) and a ordered categorical variable (e.g., an item response). Input for these later two are data frames or matrices.  Requires the mnormt package.
}
\usage{
tetrachoric(x,y=NULL,correct=.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,
     delete=TRUE)
polychoric(x,smooth=TRUE,global=TRUE,polycor=FALSE,ML=FALSE, std.err=FALSE, 
     weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,  delete=TRUE)
biserial(x,y)  
polyserial(x,y) 
polydi(p,d,taup,taud,global=TRUE,ML = FALSE, std.err = FALSE,
     weight=NULL,progress=TRUE,na.rm=TRUE,delete=TRUE,correct=.5) 
#deprecated  use polychoric instead
poly.mat(x, short = TRUE, std.err = FALSE, ML = FALSE) 
}

\arguments{
\item{x}{The input may be in one of four forms:
  
a) a data frame or matrix of dichotmous data (e.g., the lsat6 from the bock data set) or discrete numerical (i.e., not too many levels, e.g., the big 5 data set, bfi) for polychoric, or continuous for the case of biserial and polyserial. 

b) a 2 x 2 table of cell counts or cell frequencies (for tetrachoric)

c) a vector with elements corresponding to the four cell frequencies (for tetrachoric)

d) a vector with elements of the two marginal frequencies (row and column) and the comorbidity (for tetrachoric)
}
\item{y}{A (matrix or dataframe) of discrete scores. In the case of tetrachoric, these should be dichotomous, for polychoric not too many levels, for biserial they should be discrete (e.g., item responses) with not too many (<10?) categories.}
\item{correct}{Correction value to use to correct for continuity in the case of zero entry cell for tetrachoric, polychoric, polybi, and mixed.cor.  See the examples for the effect of correcting versus not correcting for continuity.}
\item{smooth}{if TRUE and if the tetrachoric matrix is not positive definite, then apply a simple smoothing algorithm using cor.smooth}

\item{global}{When finding pairwise correlations, should we use the global values of the tau parameter (which is somewhat faster), or the local values (global=FALSE)?  The local option is equivalent to the polycor solution. This will make a difference in the presence of lots of missing data. } 
\item{polycor}{A no longer used option, kept to stop other packages from breaking.}

\item{weight}{A vector of length of the number of observations that specifies the weights to apply to each case.  The NULL case is equivalent of weights of 1 for all cases.  }
\item{short}{ short=TRUE, just show the correlations, short=FALSE give the full hetcor output from John Fox's hetcor function if installed and if doing polychoric  Deprecated}
\item{std.err}{std.err=FALSE does not report the standard errors (faster)  deprecated}
\item{progress}{Show the progress bar (if  not doing multicores)}
\item{ML}{ ML=FALSE  do a quick two step procedure, ML=TRUE, do longer maximum likelihood --- very slow! Deprecated}
\item{na.rm}{Should missing data be deleted}
\item{delete}{Cases with no variance are deleted with a warning before proceeding.}
\item{p}{The polytomous input to polydi}
\item{d}{The dichotomous input to polydi}
\item{taup}{The tau values for the polytomous variables -- if global=TRUE}
\item{taud}{The tau values for the dichotomous variables -- if globabl = TRUE}

}

\details{
Tetrachoric correlations infer a latent Pearson correlation from a two x two table of frequencies with the assumption of bivariate normality.  The estimation procedure is two stage ML.  Cell frequencies for each pair of items are found. In the case of tetrachorics, cells with zero counts are replaced with .5 as a correction for continuity (correct=TRUE). 

The data typically will be a raw data matrix of responses to a questionnaire scored either true/false (tetrachoric) or with a limited number of responses (polychoric).  In both cases, the marginal frequencies are converted to normal theory thresholds and the resulting table for each item pair is converted to the (inferred)  latent  Pearson correlation that would produce the observed cell frequencies with the observed marginals.  (See \code{\link{draw.tetra}} and \code{\link{draw.cor}} for  illustrations.)

This is a very computationally intensive function which can be speeded up considerably by using multiple cores and using the parallel package.  The number of cores to use when doing polychoric or tetrachoric may be specified using the options command. The greatest step in speed is going from 1 core to 2.  This is about a 50\% savings.  Going to 4 cores seems to have about at 66\% savings, and 8 a 75\% savings.  The number of parallel processes defaults to 2 but can be modified by using the \code{\link{options}} command:  options("mc.cores"=4) will set the number of cores to 4.

The tetrachoric correlation is used in a variety of contexts, one important one being in Item Response Theory (IRT) analyses of test scores, a second in the conversion of comorbity statistics to correlation coefficients.  It is in this second context that examples of the sensitivity of the coefficient to the cell frequencies becomes apparent:

Consider the test data set from Kirk (1973) who reports the effectiveness of a ML algorithm for the tetrachoric correlation (see examples).

Examples include the lsat6 and lsat7 data sets in the \code{\link{bock}} data.

The polychoric function forms matrices of polychoric correlations by either using John Fox's polychor function or by an local function (polyc) and will also report the tau values for each alternative.  

polychoric replaces poly.mat and is recommended.   poly.mat is an alternative wrapper to the polycor function. 

biserial and polyserial correlations are the inferred latent correlations equivalent to the observed point-biserial and point-polyserial correlations (which are themselves just Pearson correlations).

The polyserial function is meant to work with matrix or dataframe input and treats missing data by finding the pairwise Pearson r corrected by the overall (all observed cases) probability of response frequency.  This is particularly useful for SAPA procedures (\url{http://sapa-project.org}) with large amounts of missing data and no complete cases.

Ability tests and personality test matrices will typically have a cleaner structure when using tetrachoric or polychoric correlations than when using the normal Pearson correlation. However, if either alpha or omega is used to find the reliability, this will be an overestimate of the squared correlation of a latent variable the observed variable.

A biserial correlation (not to be confused with the point-biserial correlation which is just a Pearson correlation) is the latent correlation between x and y where y is continuous and x is dichotomous but assumed to represent an (unobserved) continuous normal variable. Let p = probability of x level 1, and q = 1 - p.  Let zp = the normal ordinate of the z score associated with p.  Then, \eqn{rbi = r s* \sqrt(pq)/zp }.

The 'ad hoc' polyserial correlation, rps is just \eqn{r = r * sqrt(n-1)/n) \sigma y /\sum(zpi) } where zpi are the ordinates of the normal curve at the normal equivalent of the cut point boundaries between the item responses. (Olsson, 1982) 


All of these were inspired by (and adapted from) John Fox's polychor package which should be used for precise ML estimates of the correlations.  See, in particular, the hetcor function in the polychor package. 

Particularly for tetrachoric correlations from sets of data with missing data, the matrix will sometimes not be positive definite.  Various smoothing alternatives are possible, the one done here is to do an eigen value decomposition of the correlation matrix, set all negative eigen values to   10 * .Machine$double.eps, normalize the positive eigen values to sum to the number of variables, and then reconstitute the correlation matrix.   A warning is issued when this is done.

For combinations of continous, categorical, and dichotomous variables, see \code{\link{mixed.cor}}.

If using data with a variable number of response alternatives, it is necessary to use the global=FALSE option in polychoric.
}

\value{

 \item{rho}{The (matrix) of tetrachoric/polychoric/biserial correlations}
 \item{tau}{The normal equivalent of the cutpoints}

}
\references{

A. Gunther and M. Hofler. Different results on tetrachorical correlations in mplus and stata-stata announces modified procedure. Int J Methods Psychiatr Res, 15(3):157-66, 2006.

David Kirk (1973) On the numerical approximation of the bivariate normal (tetrachoric) correlation coefficient. Psychometrika, 38, 259-268.

U.Olsson, F.Drasgow, and N.Dorans (1982). The polyserial correlation coefficient. Psychometrika, 47:337-347. 

}
\author{
William Revelle
}
\note{For tetrachoric, in the degenerate case of a cell entry with zero observations, a correction for continuity is applied and .5 is added to the cell entry.  A warning is issued.  If correct=FALSE the correction is not applied.

Switched to using sadmvn from the mnormt package to speed up by 50\%. 
}



\seealso{\code{\link{mixed.cor}} to find the correlations between mixtures of continuous, polytomous, and dichtomous variables.  See also the  polychor function in the polycor package. \code{\link{irt.fa}} uses the tetrachoric function to do item analysis with the \code{\link{fa}} factor analysis function.
\code{\link{draw.tetra}} shows the logic behind a tetrachoric correlation (for teaching purpuses.)
}
\examples{
#if(require(mnormt)) {
data(bock)
tetrachoric(lsat6)
polychoric(lsat6)  #values should be the same
tetrachoric(matrix(c(44268,193,14,0),2,2))  #MPLUS reports.24

#Do not apply continuity correction -- compare with previous analysis!
tetrachoric(matrix(c(44268,193,14,0),2,2),correct=0)  

#the default is to add correct=.5 to 0 cells 
tetrachoric(matrix(c(61661,1610,85,20),2,2)) #Mplus reports .35
tetrachoric(matrix(c(62503,105,768,0),2,2)) #Mplus reports -.10
tetrachoric(matrix(c(24875,265,47,0),2,2)) #Mplus reports  0

#Do not apply continuity correction- compare with previous analysis
tetrachoric(matrix(c(24875,265,47,0),2,2), correct=0) 


#these next examples are impossible!

tetrachoric(c(0.02275000, 0.0227501320, 0.500000000))
tetrachoric(c(0.0227501320, 0.0227501320, 0.500000000))

#give a vector of two marginals and the comorbidity
tetrachoric(c(.2, .15, .1))
tetrachoric(c(.2, .1001, .1))
 #} else {
 #       message("Sorry, you must have mvtnorm installed")}

# 4 plots comparing biserial to point biserial and latent Pearson correlation
set.seed(42)
x.4 <- sim.congeneric(loads =c(.9,.6,.3,0),N=1000,short=FALSE)
y  <- x.4$latent[,1]
for(i in 1:4) {
x <- x.4$observed[,i]
r <- round(cor(x,y),1)
ylow <- y[x<= 0]
yhigh <- y[x > 0]
yc <- c(ylow,yhigh)
rpb <- round(cor((x>=0),y),2)
rbis <- round(biserial(y,(x>=0)),2)
ellipses(x,y,ylim=c(-3,3),xlim=c(-4,3),pch=21 - (x>0),
       main =paste("r = ",r,"rpb = ",rpb,"rbis =",rbis))

dlow <- density(ylow)
dhigh <- density(yhigh)
points(dlow$y*5-4,dlow$x,typ="l",lty="dashed")
lines(dhigh$y*5-4,dhigh$x,typ="l")
}




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

