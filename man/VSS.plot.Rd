\name{VSS.plot}
\alias{VSS.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot VSS fits}
\description{The Very Simple Structure criterion ( \code{\link{VSS}}) for estimating the optimal number of factors is plotted as a function of the increasing complexity and increasing number of factors.  
}
\usage{
VSS.plot(x, plottitle = "Very Simple Structure", line = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{output from VSS }
  \item{plottitle}{any title }
  \item{line}{ connect different complexities }
}
\details{Item-factor models differ in their "complexity".  Complexity 1 means that all except the greatest (absolute) loading for an item are ignored. Basically a cluster model (e.g., \code{\link{ICLUST}}). Complexity 2 implies all except the greatest two, etc.  

Different complexities can suggest different number of optimal number of factors to extract.  For personality items, complexity 1 and 2 are probably the most meaningful.

The Very Simple Structure criterion will tend to peak at the number of factors that are most interpretable for a given level of complexity.  Note that some problems, the most interpretable number of factors will differ as a function of complexity.  For instance, when doing the Harman 24 psychological variable problems, an unrotated solution of complexity one suggests one factor (g), while a complexity two solution suggests that a four factor solution is most appropriate.  This latter probably reflects a bi-factor structure.  

For examples of VSS.plot output, see \url{http://personality-project.org/r/r.vss.html}

}
\value{A plot window showing the VSS criterion varying as the number of factors and the complexity of the items.
}
\references{ \url{http://personality-project.org/r/r.vss.html}}
\author{
Maintainer: William Revelle \email{revelle@northwestern.edu}

}

\seealso{  \code{\link{VSS}},  \code{\link{ICLUST}}, \code{\link{omega}}}

\examples{
#test.data <- Harman74.cor$cov
#my.vss <- VSS(test.data)         #suggests that 4 factor complexity two solution is optimal
#VSS.plot(my.vss)                 #see the graphics window

## The function is currently defined as
function(x,plottitle="Very Simple Structure",line=FALSE)
  {
  	n=dim(x)
  	symb=c(49,50,51,52)              #plotting sym
	plot(x$cfit.1,ylim=c(0,1),type="b",ylab="Very Simple Structure Fit",xlab="Number of Factors",pch=49)
    if (line) lines(x$fit)
	 

	title(main=plottitle)
	x$cfit.2[1]<-NA
	x$cfit.3[1]<-NA
	x$cfit.3[2]<-NA
	x$cfit.4[1]<-NA
	x$cfit.4[2]<-NA
	x$cfit.4[3]<-NA
	lines(x$cfit.2)
	points(x$cfit.2,pch=50)
	lines(x$cfit.3)
	points(x$cfit.3,pch=symb[3])
	lines(x$cfit.4)
	points(x$cfit.4,pch=symb[4])
  }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
