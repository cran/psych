\name{mat.regress}
\alias{mat.regress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multiple Regression from matrix input}
\description{
  This function  extracts subsets of variables (x and y) from a correlation 
  matrix (m) and then find the multiple correlation and beta weights of the (x) set 
  predicting each member of the (y) set.  
}
\usage{
mat.regress(m, x, y,digits=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{  a matrix of correlations }
  \item{x}{  the column numbers of the x set (e.g., c(1,3,5)}
  \item{y}{  the column numbers of the y set  (e.g., c(2,4,6)}
  \item{digits}{round the answer to digits}
}
\details{
  Although it is more common to calculate multiple regression from raw data, it is, 
  of course, possible to do so from a set of correlations.  The input to the function is 
  a square covariance or correlation matrix, as well as the column numbers of the x (predictor) and 
  y (criterion) variables. 
  
  The output is a set of multiple correlations, one for each dependent variable in the y set. 
  
  A typical use in the SAPA project is to form item composites by clustering or factoring (see \code{\link{ICLUST}}, \code{\link{principal}}), extract the clusters from these results (\code{\link{factor2cluster}}), and then form the composite correlation matrix using \code{\link{cluster.cor}}.  The variables in this reduced matrix may then be used in multiple R procedures using mat.regress.
}
\value{
  
  \item{beta }{the beta weights for each variable in X for each variable in Y}
  \item{R2 }{The multiple R2 for each equation}
 
}
\references{ For an application of this procedure, see http://personality-project.org/revelle/publications/sapa.pdf }
\author{William Revelle \cr
        Department of Psychology \cr
        Northwestern University  \cr
        Evanston, Illiniois  \cr


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
#VisualPerception 0.397               0.22    0.162        0.296
#Cubes            0.064               0.18    0.056        0.049
#PaperFormBoard   0.125               0.10   -0.158        0.005
#
#$R2
#             Flags GeneralInformation           Addition       CountingDots 
#             0.239              0.148              0.034              0.101 
#


## The function is currently defined as
function(m,x,y)  {
 #a function to extract subsets of variables (a and b) from a correlation matrix m
  #and find the multiple correlation beta weights + R2 of the a set predicting the b set
    
    #first reorder the matrix to select the right variables
         nm <- dim(m)[1]
        t.mat <- matrix(0,ncol=nm,nrow=nm)
        xy <- c(x,y)
         numx <- length(x)
     	numy <- length(y)
        nxy <- numx+numy
        for (i in 1:nxy) {
     	t.mat[i,xy[i]] <- 1 }
     	
     	reorder <- t.mat \%*\% m \%*\% t(t.mat)
     	a.matrix <- reorder[1:numx,1:numx]
     	b.matrix <- reorder[1:numx,(numx+1):nxy]
     	model.mat <- solve(a.matrix,b.matrix)       #solve the equation bY~aX
     	if (length(y) >1 ) { rownames(model.mat) <- rownames(m)[x]
     	 colnames(model.mat) <- colnames(m)[y]
     	 
     	R2 <- colSums(model.mat * b.matrix) }
     	 else { R2 <- sum(model.mat * b.matrix)
     	 names(model.mat) <- rownames(m)[x]
     	 names(R2) <- colnames(m)[y]}
     	mat.regress <- list(beta=model.mat,R2=R2)
     	return(mat.regress)
     	}
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line
