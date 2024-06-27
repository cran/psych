\name{RV}
\alias{RV}
\title{Three measures of the correlations  between sets of variables}

\description{How to measure the the correlation between two clusters or groups of variables x and y from the same data set is a recurring problem.  Perhaps the most obvious is simply the unweighted correlation Ru. 

Consider the matrix M composed of four submatrices

 \tabular{lll}{

    \tab   Rx  \tab Rxy \cr
 M = \tab  Rxy  \tab Ry  \cr
 }  



The unit weighted correlation, Ru is merely 
 
  \eqn{Ru =\frac{\Sigma{r_{xy}}}{\sqrt{\Sigma{r_x}\Sigma{r_y}} }}{Ru = \Sigma{r_{xy}/\sqrt{\Sigma{r_x}}\Sigma{r_y}}}
  \cr 
  
 A second is the Set correlation (also found in \code{\link{lmCor}}) by Cohen 1982) which is  \cr
\eqn{Rset = 1- \frac{det(m)}{det(x)* det(y)}}
 \cr
 Where m is the full matrix (x+y)by (x+y). and det represents the determinant.  
 \cr
 
 A third approach (the RV coeffiecent) was introduced by  Escoufier (1970) and Robert and Escoufier (1976).
 \cr
 \eqn{RV = \frac{tr(xy (xy)')}{\sqrt{(tr(x x') * tr(y  
 y'))}}}{RV=   tr( xy  \%*\% t(xy))/\sqrt{(tr(x \%*\% t(x)))* (tr(y \%*\% t(y))}}.
   
 Where \code{\link{tr}} is the trace operator.  (The sum of the diagonals).
 
  The analysis can be done from the raw data or from correlation or covariance matrices. From the raw data, just specify the x and y variables.  If using correlation/covariance matrixes,  the xy matrix must be specified as well.   
 
 If using raw data, just specify the x and y columns and the data file.  
}

\usage{
RV(x, y, xy = NULL, data=NULL, cor = "cor",correct=0)
}
\arguments{
  \item{x}{ Columns of the data matrix of n rows and p columns, (if data is specified) or a p x p correlation matrix.
}
  \item{y}{Columns of a raw data matrix of n rows and q columns, or a q * q correlation matrix.

}
  \item{xy}{A p x q correlation or covariance matrix, if not using the raw data.
}
\item{data}{A matrix or data frame containing the raw data. }

  \item{cor}{If xy is NULL,  find the p x p correlations or covariances from x, and the q x q correlations from  y as  well as the p x q covariance/correlation matrix..  Options are "cor" (for Pearson), "spearman" , "cov" for covariances, "tet" for tetrachoric, or "poly" for polychoric correlation.}
  
  \item{correct}{The correction for continuity if desired.  }
}

\details{If using raw data, just specify the  columns in x and y.  If using a correlation matrix or covariance matrix, the inter corrlations/covariances) are specified in xy.  The results match those of the RV function from matrixCalculations and the coeffRV function from  factoMineR. 
}
 
\value{
\item{RV}{The RV statistic}
\item{Rset}{Cohen's set correlation}
\item{Ru}{The unit weighted correlation between x and y.  }
\item{Rx}{The correlation matrix of the x variables.}
\item{Ry}{The correlation matrix of the y variables.}
\item{Rxy}{The intercorrelations of x and y.}
}


\references{
P. Robert and Y. Escoufier, 1976, A Unifying Tool for Linear Multivariate Statistical Methods: The RV- Coefficient. Journal of the Royal Statistical Society. Series C (Applied Statistics), Volume 25, pp. 257-265.

J. Cohen (1982)  Set correlation as a general multivariate data-analytic method. Multivariate Behavioral Research, 17(3):301-341.

}
\author{William Revelle}

\seealso{\code{\link{lmCor}}  for unit weighted correlations. 

}
\examples{
#from raw data
RV (attitude[1:3],attitude[4:7])  #find the correlations
RV (attitude[1:3],attitude[4:7],cor="cov")

#find the correlations
R <- cor(attitude)
r1 <- R[1:3,1:3]
r2 <- R[4:7,4:7]
r12 <- R[1:3,4:7]
RV(r1,r2,r12)

#or find the covariances
C <- cov(attitude)
c1 <- C[1:3,1:3]
c2 <- C[4:7,4:7]
c12 <- C[1:3,4:7]
RV(c1, c2, c12)

}
\keyword{ multivariate }
\keyword{models}