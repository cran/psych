
\name{congruence}
\alias{congruence}
\alias{cohen.profile}
\alias{distance}
\title{Matrix and profile congruences and distances}
\description{The congruence coefficient two matrices is just the cross product of their respective values divided by the square root of their sums of squares. If the columns are zero centered, this is just the correlation. If the columns are centered around the scale neutral point, this is Cohen's profile correlation. A set of distances (city block, euclidean, Minkowski) may be found by the distance function. 
}
\usage{
congruence(x,y=NULL)
cohen.profile(x,y=NULL ,M=NULL)
distance(x,y=NULL,r=2)
}
\arguments{
   \item{x}{ A matrix of factor loadings or a list of matrices of factor loadings}
  \item{y}{ A second matrix of factor loadings (if x is a list, then y may be empty)}
  
\item{M}{The midpoint of the items which should be used for reflection. NULL if to be found from the data. }
  \item{r}{The exponent for the generalized distance. (r=1 is city block, r=2 is euclidian, r=100 or larger emphasize the largest distance )}
}
\details{
Congruences are the cosines of pairs of vectors defined by a matrix and based at the origin.  Thus, for values that differ only by a scaler the congruence will be 1.

For two matrices, F1 and F2, the measure of  congruence, phi, is 
\deqn{
\phi = \frac{\sum F_1 F_2}{\sqrt{\sum(F_1^2)\sum(F_2^2)}} 
.}{phi = sum(F1 F2)/sqrt(sum(F1^2) sum(F2^2)) }

It is an interesting exercise to compare  congruences with the correlations of the respective values. Congruences are based upon the raw cross products, while correlations are based upon centered cross products. That is,  correlations of items are cosines of the vectors based at the mean loading for each column.   

\deqn{
\phi = \frac{\sum (F_1-a) (F_2 - b)}{\sqrt{\sum((F_1-a)^2)\sum((F_2-b)^2)}} 
.}{phi = sum((F1-a)(F2-b))/sqrt(sum((F1-a)^2) sum((F2-b)^2)) }.

For congruence coefficients, a = b= 0.  For correlations a=mean F1, b= mean F2.

Input may either be one or two matrices.

For congruences of factor or component loading matrices, use \code{\link{factor.congruence}}.

Normally, all factor loading matrices should be complete (have no missing loadings).  In the case where some loadings are missing, if the use option is specified, then variables with missing loadings are dropped.


If the data are zero centered, this is the correlation, if the data are centered around the scale midpoint (M), this is Cohen's Similarity coefficient.  See examples. If M is not specified, it is found as the midpoint of the items in x and y.

\code{\link{cohen.profile}} applies the \code{\link{congruence}} function to data centered around M.  M may be specified, or found from the data. The last example is taken from Cohen (1969). 

\code{\link{distance}} finds the generalized distance as a function of r.  City block (r=1), Euclidean (r=2) or weighted towards maximimum (r >2).

}
\value{A matrix of  congruences or distances.
}

\references{
Burt, Cyril (1948) Methods of factor-analysis with and without successive approximation.  British Journal of Educational Psychology, 7 (2) 172-195.

Burt, Cyril (1948) The factorial study of temperamental traits. British Journal of Statistical Psychology, 1(3) 178-203.

Cohen, Jacob (1969), rc: A profile similarity coefficient invariant over variable reflection. Psychological Bulletin, 71 (4) 281-284.

}
\author{William Revelle
}


\seealso{\code{\link{factor.congruence}}.

}
\examples{#cohen's example
# a and b have reversed one item around the midpoint
co <-  data.frame(ira=c(2,6,5,6,4),
       jim=c(1,3,5,4,4),
      a=c(5,6,5,6,4),b=c(6,3,5,4,4))

lowerMat(congruence(co-3.5)) # congruence around the midpoint is insensitive to reflection
lowerCor(co)   #the correlation is not
lowerMat(congruence(scale(co,scale=FALSE))) #zero centered congruence is  r
cohen.profile(co)
}
\keyword{ multivariate }
\keyword{ models }

% Use only one concept per line.
