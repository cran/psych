\name{Yule}
\alias{Yule}
\alias{Yule.inv}
\alias{Yule2phi}
\alias{Yule2tetra}
\alias{Yule2poly}
\title{From a two by two table, find the Yule coefficients of association, convert to phi, or tetrachoric, recreate table the table to create the Yule coefficient.}
\description{One of the many measures of association is the Yule coefficient.  Given a two x two table of counts \cr
\tabular{llll}{
\tab a \tab b \tab R1 \cr
\tab c \tab d \tab R2 \cr
\tab C1 \tab C2 \ tab n \cr
}
Yule Q is (ad - bc)/(ad+bc). \cr
Conceptually, this is the number of pairs in agreement (ad) - the number in disagreement (bc) over the total number of paired observations.  Warren (2008) has shown that  Yule's Q is one of the ``coefficients that have zero value under statistical independence, maximum value unity, and minimum value minus unity independent of the marginal distributions" (p 787). 
\cr
ad/bc is the odds ratio and Q = (OR-1)/(OR+1) 
\cr
Yule's coefficient of colligation is Y = (sqrt(OR) - 1)/(sqrt(OR)+1)
Yule.inv finds the cell entries for a particular Q and the marginals (a+b,c+d,a+c, b+d).  This is useful for converting old tables of correlations into more conventional \code{\link{phi}} or tetrachoric correlations \code{\link{tetrachoric}}
\cr
Yule2phi and Yule2tetra convert the Yule Q with set marginals to the correponding phi or tetrachoric correlation.
}
\usage{
Yule(x,Y=FALSE)  #find Yule given a two by two table of frequencies
 #find the frequencies that produce a Yule Q given the Q and marginals
Yule.inv(Q,m,n=NULL)   
#find the phi coefficient that matches the Yule Q given the marginals
Yule2phi(Q,m,n=NULL)    
Yule2tetra(Q,m,n=NULL,correct=TRUE) 
   #Find the tetrachoric correlation given the Yule Q and the marginals
#(deprecated) Find the tetrachoric correlation given the Yule Q and the marginals   
Yule2poly(Q,m,n=NULL,correct=TRUE)   
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of four elements or a two by two matrix }
  \item{Y}{Y=TRUE return Yule's Y coefficient of colligation}
  \item{Q}{Either a single Yule coefficient or a matrix of Yule coefficients}
  \item{m}{The vector c(R1,C2) or a two x two matrix of marginals or a four element vector of marginals.  The preferred form is c(R1,C1)}
  \item{n}{The number of subjects (if the marginals are given as frequencies}
  \item{correct}{When finding a tetrachoric correlation, should small cell sizes be corrected for continuity.  See \code{\{link{tetrachoric}} for a discussion.}
}
\details{Yule developed two measures of association for two by two tables.  Both are functions of the odds ratio 
}
\value{
  \item{Q}{The Yule Q coefficient}
  \item{R}{A two by two matrix of counts}
  \item{result}{If given matrix input, then a matrix of phis or tetrachorics}
}
\references{Yule, G. Uday (1912) On the methods of measuring association between two attributes. Journal of the Royal Statistical Society, LXXV, 579-652

Warrens, Matthijs (2008), On Association Coefficients for 2x2 Tables and Properties That Do Not Depend on the Marginal Distributions. Psychometrika, 73, 777-789. 

}

\author{ William Revelle }
\note{Yule.inv is currently done by using the optimize function, but presumably could be redone by solving a quadratic equation.
}
\seealso{ See Also as \code{\link{phi}}, \code{\link{tetrachoric}},  \code{\link{Yule2poly.matrix}}, \code{\link{Yule2phi.matrix}} }
\examples{
Nach <- matrix(c(40,10,20,50),ncol=2,byrow=TRUE)
Yule(Nach)
Yule.inv(.81818,c(50,60),n=120)
Yule2phi(.81818,c(50,60),n=120)
Yule2tetra(.81818,c(50,60),n=120)
phi(Nach)  #much less
#or express as percents and do not specify n
Nach <- matrix(c(40,10,20,50),ncol=2,byrow=TRUE)
Nach/120
Yule(Nach)
Yule.inv(.81818,c(.41667,.5))
Yule2phi(.81818,c(.41667,.5))
Yule2tetra(.81818,c(.41667,.5))
phi(Nach)  #much less
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{models }% __ONLY ONE__ keyword per line
