\name{mixed.cor}
\alias{mixed.cor}
\title{Find correlations for mixtures of continuous, polytomous, and dichotomous variables}
\description{For data sets with continuous, polytomous and dichotmous variables, the absolute Pearson correlation is downward biased from the underlying latent correlation.  mixed.cor finds Pearson correlations for the continous variables, \code{\link{polychoric}}s for the polytomous items, \code{\link{tetrachoric}}s for the dichotomous items, and the \code{\link{polyserial}} or \code{\link{biserial}} correlations for the various mixed variables. Results include the complete correlation matrix, as well as the separate correlation matrices and difficulties for the polychoric and tetrachoric correlations.
}
\usage{
mixed.cor(x = NULL, p = NULL, d = NULL,smooth=TRUE,ncat=8)
}
\arguments{
  \item{x}{A set of continuous variables (may be missing) or, if p and d are missing, the variables to be analyzed.}
  \item{p}{A set of polytomous items (may be missing)}
  \item{d}{A set of dichotomous items (may be missing)}
  \item{smooth}{If TRUE, then smooth the correlation matix if it is non-positive definite}
  \item{ncat}{The number of categories beyond which a variable is considered "continuous".}
}
\details{
This function is particularly useful as part of the Synthetic Apeture Personality Assessment (SAPA) data sets where continuous variables (age, SAT V, SAT Q, etc) and mixed with polytomous personality items taken from the International Personality Item Pool (IPIP) and the dichotomous experimental IQ items that have been developed as part of SAPA (see, e.g., Revelle, Wilt and Rosenthal, 2010).  

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
\note{Note that the variables must be organized by type of data: first continuous, then polytomous, then dichotomous.  }

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
lower.mat(rp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}