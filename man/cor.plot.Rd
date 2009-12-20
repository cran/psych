\name{cor.plot}
\Rdversion{1.1}
\alias{cor.plot}
\title{Create an image plot for a correlation or factor matrix}

\description{Correlation matrices may be shown graphically by using the image function to emphasize structure.  This is a particularly useful tool for showing the structure of small correlation matrices with a clear structure.  Meant for the pedagogical value of the graphic for teaching or discussing factor analysis and other multivariate techniques.
}
\usage{
cor.plot(r,colors=FALSE, n=10,main=NULL,zlim=c(0,1),show.legend=TRUE,...)
}
\arguments{
  \item{r}{A correlation matrix or the output of  \code{\link{factor.pa}}, \code{\link{factor.minres}} or \code{\link{omega}}.  }
  \item{colors}{Defaults to FALSE (grey), but colors=TRUE will use topo.colors}
  \item{n}{The number of levels of shading to use.  Defaults to 10}
  \item{main}{A title. Defaults to ``correlation plot"}
  \item{zlim}{The range of values to color -- defaults to 0 to 1}
  \item{show.legend}{A legend (key) to the colors is shown on the right hand side}
  \item{...}{Other parameters for axis (e.g., cex.axis to change the font size)}
}
\details{When teaching about factor analysis or cluster analysis, it is useful to graphically display the structure of correlation matrices.  This is a simple graphical display using the image function. 

The difference of mat.plot with a regular image plot is that the primary diagonal goes from the top left to the lower right. 

The zlim parameter defaults to 0 to 1.  This means that negative correlations are treated as zero.  This is advantageous when showing general factor structures, because it makes the 0 white.  

The default shows a legend for the color coding on the right hand side of the figure.

Inspired, in part, by a paper by S. Dray (2008)  on the number of components problem. 
}

\author{William Revelle
}

\references{
Dray, Stephane (2008) On the number of principal components: A test of dimensionality based on measurements of similarity between matrices. Computational Statistics \& Data Analysis. 52, 4, 2228-2237.
}

\examples{
data(bifactor)
cor.plot(Thurstone,TRUE, main="9 cognitive variables from Thurstone")
simp <- sim.circ(24)
cor.plot(cor(simp),colors=TRUE,zlim=c(-1,1),main="24 variables in a circumplex")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line

