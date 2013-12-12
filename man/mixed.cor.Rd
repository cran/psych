\name{mixed.cor}
\alias{mixed.cor}
\title{Find correlations for mixtures of continuous, polytomous, and dichotomous variables}
\description{For data sets with continuous, polytomous and dichotmous variables, the absolute Pearson correlation is downward biased from the underlying latent correlation.  mixed.cor finds Pearson correlations for the continous variables, \code{\link{polychoric}}s for the polytomous items, \code{\link{tetrachoric}}s for the dichotomous items, and the \code{\link{polyserial}} or \code{\link{biserial}} correlations for the various mixed variables. Results include the complete correlation matrix, as well as the separate correlation matrices and difficulties for the polychoric and tetrachoric correlations.
}
\usage{
mixed.cor(x = NULL, p = NULL, d=NULL,smooth=TRUE, correct=TRUE,global=TRUE, 
        ncat=8, polycor=FALSE,use="pairwise",method="pearson",weight=NULL)
}
\arguments{
  \item{x}{A set of continuous variables (may be missing) or, if p and d are missing, the variables to be analyzed.}
  \item{p}{A set of polytomous items (may be missing)}
  \item{d}{A set of dichotomous items (may be missing)}
  \item{smooth}{If TRUE, then smooth the correlation matix if it is non-positive definite}
  \item{correct}{When finding tetrachoric correlations, should 0 values be corrected for continuity?}
  \item{global}{For polychorics, should the global values of the tau parameters be used, or should the pairwise values be used.  Set to local if errors are occurring.}
  \item{ncat}{The number of categories beyond which a variable is considered "continuous".}
  \item{polycor}{Should the polycor function from the polychor package be used?}
  \item{use}{The various options to the \code{\link{cor}} function include "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". The default here is "pairwise"}
  \item{method}{The correlation method to use for the continuous variables. "pearson" (default), "kendall", or "spearman"}
  \item{weight}{If specified, this is a vector of weights (one per participant) to differentially weight participants. The NULL case is equivalent of weights of 1 for all cases. }
}
\details{
This function is particularly useful as part of the Synthetic Apeture Personality Assessment (SAPA) (\url{http://sapa-project.org}) data sets where continuous variables (age, SAT V, SAT Q, etc) and mixed with polytomous personality items taken from the International Personality Item Pool (IPIP) and the dichotomous experimental IQ items that have been developed as part of SAPA (see, e.g., Revelle, Wilt and Rosenthal, 2010).  

Item response analyses using \code{\link{irt.fa}} may be done separately on the polytomous and dichotomous items  in order to develop internally consistent scales. These scale may, in turn, be correlated with each other using the complete correlation matrix found by mixed.cor and using the \code{\link{score.items}} function.

This function is not quite as flexible as the hetcor function in John Fox's polychor package. 

Note that the variables may be organized by type of data: first continuous, then polytomous, then dichotomous.  This is advantageous  }
\value{
 \item{rho}{The complete matrix}
 \item{rx}{The Pearson correlation matrix for the continuous items}
 \item{poly}{the polychoric correlation (poly$rho) and the item difficulties (poly$tau)}
 \item{tetra}{the tetrachoric correlation (tetra$rho) and the item difficulties (tetra$tau)}
}

\author{William Revelle}
\note{mixed.cor was designed for the SAPA project (\url{http://sapa-project.org}) with large data sets with a mixture of continuous, dichotomous, and polytomous data.  For smaller data sets, it is sometimes the case that the global estimate of the tau parameter will lead to unstable solutions.  This may be corrected by setting the global parameter = FALSE.  

When finding correlations between dummy coded SAPA data (e.g., of occupations), the real correlations are all slightly less than zero because of the ipsatized nature of the data.  This leads to a non-positive definite correlation matrix because the matrix is no longer of full rank. Smoothing will correct this, even though this might not be desired.  Turn off smoothing in this case.

 Note that the variables no longer need to be organized by type of data: first continuous, then polytomous, then dichotomous.  }

\references{
W.Revelle, J.Wilt, and A.Rosenthal. Personality and cognition: The personality-cognition link. In A.Gruszka, G. Matthews, and B. Szymura, editors,  Handbook of Individual Differences in Cognition: Attention, Memory and Executive
  Control, chapter 2, pages 27-49. Springer, 2010.}



\seealso{ \code{\link{polychoric}}, \code{\link{tetrachoric}}, \code{\link{score.items}}, \code{\link{score.irt}}
}
\examples{
data(bfi) 
r <- mixed.cor(bfi[,c(1:5,26,28)])
r
#compare to raw Pearson
#note that the biserials and polychorics are not attenuated
rp <- cor(bfi[c(1:5,26,28)],use="pairwise")
lowerMat(rp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}