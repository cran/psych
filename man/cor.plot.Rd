\name{cor.plot}
\Rdversion{1.1}
\alias{cor.plot}
\title{Create an image plot for a correlation or factor matrix}

\description{Correlation matrices may be shown graphically by using the image function to emphasize structure.  This is a particularly useful tool for showing the structure of  correlation matrices with a clear structure.  Partially meant for the pedagogical value of the graphic for teaching or discussing factor analysis and other multivariate techniques.
}
\usage{
cor.plot(r,numbers=FALSE,colors=TRUE, n=51,main=NULL,zlim=c(-1,1),show.legend=TRUE,labels=NULL,n.legend=10,keep.par=TRUE,...)
}
\arguments{
  \item{r}{A correlation matrix or the output of  \code{\link{factor.pa}}, \code{\link{factor.minres}} or \code{\link{omega}}.  }
  \item{numbers}{Display the numeric value of the correlations.  Defaults to FALSE.}
  \item{colors}{Defaults to TRUE and colors use colors from the colorRampPalette from red through white to blue, but colors=FALSE will use a grey scale}
  \item{n}{The number of levels of shading to use.  Defaults to 51}
  \item{main}{A title. Defaults to "correlation plot"}
  \item{zlim}{The range of values to color -- defaults to -1 to 1}
  \item{show.legend}{A legend (key) to the colors is shown on the right hand side}
  \item{labels}{if NULL, use column and row names, otherwise use labels}
  \item{n.legend}{How many categories should be labelled in the legend?}
  \item{keep.par}{restore the graphic parameters when exiting}
  \item{...}{Other parameters for axis (e.g., cex.axis to change the font size)}
}
\details{When summarizing the correlations of large data bases or when teaching about factor analysis or cluster analysis, it is useful to graphically display the structure of correlation matrices.  This is a simple graphical display using the image function. 

The difference between mat.plot with a regular image plot is that the primary diagonal goes from the top left to the lower right. 
zlim defines how to treat the range of possible values. -1 to 1 and the color choice is more reasonable.  Setting it as c(0,1) will lead to negative correlations  treated as zero.  This is advantageous when showing general factor structures, because it makes the 0 white.  

The default shows a legend for the color coding on the right hand side of the figure.

Inspired, in part, by a paper by S. Dray (2008)  on the number of components problem. 

Modified following suggestions by David Condon and Josh Wilt to use a more meaningful color choice ranging from dark red (-1) through white (0) to dark blue (1).  Further modified to include the numerical value of the correlation.  (Inspired by the corrplot package).
}

\author{William Revelle
}

\references{
Dray, Stephane (2008) On the number of principal components: A test of dimensionality based on measurements of similarity between matrices. Computational Statistics \& Data Analysis. 52, 4, 2228-2237.
}

\seealso{ \code{\link{fa}}, \code{\link{mat.sort}}
}
\examples{
cor.plot(Thurstone,main="9 cognitive variables from Thurstone") #just blue implies positive manifold
cor.plot(Thurstone, zlim=c(0,1),main="9 cognitive variables from Thurstone")
cor.plot(mat.sort(Thurstone),TRUE,zlim=c(0,1), main="9 cognitive variables from Thurstone (sorted by factor loading) ")
simp <- sim.circ(24)
cor.plot(cor(simp),main="24 variables in a circumplex")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line

