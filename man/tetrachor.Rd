\name{tetra}
\alias{tetrachor}
\alias{tetra}

\title{Tetrachoric correlations from various types of input}

\description{The tetrachoric correlation is the estimate of a Pearson Correlation from a two x two table with the assumption of bivariate normality. Particularly important when doing Item Response Theory or converting comorbidity statistics normal theory correlations. Input may be a 2 x 2 table of cell frequencies, a vector of cell frequencies, or a data.frame or matrix of dichotomous data.
}
\usage{
tetrachor(x,correct=TRUE)
}

\arguments{
  \item{x}{The input may be in one of four forms:
  
a)a data frame or matrix of dichotmous data (e.g., the lsat6 from the bock data set)

b)a 2 x 2 table of cell counts or cell frequencies

c) a vector with elements corresponding to the four cell frequencies

d) a vector with elements of the two marginal frequencies (row and column) and the comorbidity
} 
  \item{correct}{Correct for continuity in the case of zero entry cell}
}

\details{
Tetrachoric correlations convert a two x two table of frequencies to a Pearson correlation with the assumption of bivariate normality.  The estimation procedure is two stage ML.  Cells with zero counts are replaced with .5 as a correction for continuity.

The tetrachoric correlation is used in a variety of contexts, one important one being in Item Response Theory (IRT) analyses of test scores, a second in the conversion of comorbity statistics to correlation coefficients.  It is in this second context that examples of the sensitivity of the coefficient to the cell frequencies becomes apparent:

Consider the test data set from Kirk (1973) who reports the effectiveness of a ML algorithm for the tetrachoric correlation (see examples).

Examples include the lsat6 and lsat7 data sets in the bock data.
}

\value{

 \item{rho}{The (matrix) of tetrachoric correlations}
 \item{tau}{The normal equivalent of the cutpoints}

}
\references{

A. Gunther and M. Hofler. Different results on tetrachorical correlations in mplus and stata-stata announces modified procedure. Int J Methods Psychiatr Res, 15(3):157-66, 2006.

David Kirk (1973) On the numerical approximation of the bivariate normal (tetrachoric) correlation coefficient. Psychometrika, 38, 259-268.

}
\author{
William Revelle
}
\note{In the degenerate case of a cell entry with zero observations, a correction for continuity is applied and .5 is added to the cell entry.  A warning is issued.  
}



\seealso{ See also the  polychor function in the polycor package.
}
\examples{
if(require(mvtnorm)) {
data(bock)
tetrachor(lsat6)
tetrachor(matrix(c(44268,193,14,0),2,2))  #MPLUS reports.24
tetrachor(matrix(c(61661,1610,85,20),2,2)) #Mplus reports .35
tetrachor(matrix(c(62503,105,768,0),2,2)) #Mplus reports -.10
tetrachor(matrix(c(24875,265,47,0),2,2)) #Mplus reports  0
tetrachor(c(0.02275000, 0.0227501320, 0.500000000))
tetrachor(c(0.0227501320, 0.0227501320, 0.500000000)) } else {message("Sorry, you must have mvtnorm installed")}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

