\name{mat.regress}
\alias{mat.regress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multiple Regression from matrix input}
\description{
 Extract subsets of variables (x and y) from a correlation or data matrix
  matrix and find the multiple correlation and beta weights of the (x) set 
  predicting each member of the (y) set.  
}
\usage{
mat.regress(m, x, y,n.obs=NULL,digits=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{  a matrix of correlations or, if not square of data}
  \item{x}{  the column numbers of the x set (e.g., c(1,3,5)}
  \item{y}{  the column numbers of the y set  (e.g., c(2,4,6)}
  \item{n.obs}{If specified, then confidence intervals, etc. are calculated, not needed if raw data are given}
  \item{digits}{round the answer to digits}
}
\details{
  Although it is more common to calculate multiple regression from raw data, it is, 
  of course, possible to do so from a set of correlations.  The input to the function is 
  a square covariance or correlation matrix, as well as the column numbers of the x (predictor) and 
  y (criterion) variables. The program will find correlations if given raw data.
  
  The output is a set of multiple correlations, one for each dependent variable in the y set. 
  
  A typical use in the SAPA project is to form item composites by clustering or factoring (see \code{\link{ICLUST}}, \code{\link{principal}}), extract the clusters from these results (\code{\link{factor2cluster}}), and then form the composite correlation matrix using \code{\link{cluster.cor}}.  The variables in this reduced matrix may then be used in multiple R procedures using mat.regress.

 Although the overall matrix can have missing correlations, the correlations in the subset of the matrix used for prediction must exist.
 
 If the number of observations is entered, then the conventional confidence intervals, statistical significance, and shrinkage estimates are  reported.
 
 If the input matrix is rectangular, correlations are found from the data.
}
\value{
  
  \item{beta }{the beta weights for each variable in X for each variable in Y}
  \item{R}{The multiple R for each equation (the amount of change a unit in the predictor set leads to in the criterion set). }
  \item{R2 }{The multiple R2 (\% variance acounted for) for each equation}
 
}
\author{William Revelle \cr
       
Maintainer: William Revelle <revelle@northwestern.edu> }


\seealso{\code{\link{cluster.cor}}, \code{\link{factor2cluster}},\code{\link{principal}},\code{\link{ICLUST}} }
\examples{

\dontrun{
test.data <- Harman74.cor$cov     #24 mental variables
#choose 3 of them to regress against another 4 -- arbitrary choice of variables
print(mat.regress(test.data,c(1,2,3),c(4,5,10,12)),digits=2)
}
#gives this output
#print(mat.regress(test.data,c(1,2,3),c(4,5,10,12)),digits=2)
#$beta
#                 Flags GeneralInformation Addition CountingDots
#VisualPerception  0.40               0.22     0.16         0.30
#Cubes             0.06               0.18     0.06         0.05
#PaperFormBoard    0.12               0.10    -0.16         0.00
#
#$R
#             Flags GeneralInformation           Addition       CountingDots 
#             0.49               0.38               0.18               0.32 
#
#$R2
#             Flags GeneralInformation           Addition       CountingDots 
#              0.24               0.15               0.03               0.10 
#
#
data(attitude)
mat.regress(attitude,c(1:3),c(4:7)) #standardized regression from raw data
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line
