\name{omega.graph}
\alias{omega.diagram}
\alias{omega.graph}
\title{Graph hierarchical factor structures }
\description{Hierarchical factor structures represent the correlations between variables in terms of a smaller set of correlated factors which themselves can be represented by a higher order factor.

Two alternative solutions to such structures are found by the \code{\link{omega}} function.  The correlated factors solutions represents the effect of the higher level, general factor, through its effect on the correlated factors.  The other representation makes use of the Schmid Leiman transformation to find the direct effect of the general factor upon the original variables as well as the effect of orthogonal residual group factors upon the items. 

Graphic presentations of these two alternatives are helpful in understanding the structure.  omega.graph  and omega.diagram draw both such structures.  Graphs are drawn directly onto the graphics window or expressed in ``dot" commands for conversion to graphics using implementations of Graphviz (if using omega.graph).

Using Graphviz allows the user to clean up the Rgraphviz output. However, if Graphviz and Rgraphviz are not available, use omega.diagram.  

See the other structural diagramming functions, \code{\link{fa.diagram}} and \code{\link{structure.diagram}}. 

In addition 
}
\usage{
omega.diagram(om.results,sl=TRUE,sort=TRUE,labels=NULL,flabels=NULL,cut=.2,
gcut=.2,simple=TRUE,  errors=FALSE, digits=1,e.size=.1,rsize=.15,side=3,
    main=NULL,cex=NULL,color.lines=TRUE,marg=c(.5,.5,1.5,.5),adj=2, ...) 
omega.graph(om.results, out.file = NULL,  sl = TRUE, labels = NULL, size = c(8, 6), 
    node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10),  
    rank.direction=c("RL","TB","LR","BT"), digits = 1, title = "Omega", ...)
}

\arguments{
  \item{om.results}{The output from the omega function }
  \item{out.file}{ Optional output file for off line analysis using Graphviz }
  \item{sl}{ Orthogonal clusters using the Schmid-Leiman transform (sl=TRUE) or oblique clusters }
  \item{labels}{ variable labels }
  \item{flabels}{Labels for the factors (not counting g)}
  \item{size}{size of graphics window }
  \item{node.font}{ What font to use for the items}
  \item{edge.font}{What font to use for the edge labels }
  \item{rank.direction}{ Defaults to left to right }
  \item{digits}{ Precision of labels }
  \item{cex}{control font size}
  \item{color.lines}{Use black for positive, red for negative}
  \item{marg}{The margins for the figure are set to be wider than normal by default}
  \item{adj}{Adjust the location of the factor loadings to vary as factor mod 4 + 1}
  \item{title}{ Figure title }
    \item{main}{ main figure caption }
  \item{\dots}{Other options to pass into the graphics packages }
  \item{e.size}{the size to draw the ellipses for the factors. This is scaled by the  number of variables.}
  \item{cut}{Minimum path coefficient to draw}
  \item{gcut}{Minimum general factor path to draw}
  \item{simple}{draw just one path per item}
  \item{sort}{sort the solution before making the diagram}
  \item{side}{on which side should errors  be drawn?}
  \item{errors}{show the error estimates}
  \item{rsize}{size of the rectangles}
}
\details{While omega.graph requires the Rgraphviz package, omega.diagram does not.  \code{\link{omega}} requires the GPArotation package.

 
}
\value{
  
  \item{clust.graph }{A graph object}
  \item{sem}{A matrix suitable to be run throughe the sem function in the sem package.}
 
 }
\references{ \url{https://personality-project.org/r/r.omega.html} \cr

Revelle, W. (in preparation) An Introduction to  Psychometric Theory with applications in R.  \url{https://personality-project.org/r/book/}

Revelle, W. (1979).  Hierarchical cluster analysis and the internal structure of tests. Multivariate Behavioral Research, 14, 57-74. (\url{https://personality-project.org/revelle/publications/iclust.pdf})

 Zinbarg, R.E., Revelle, W., Yovel, I., & Li. W.  (2005). Cronbach's Alpha, Revelle's Beta, McDonald's Omega: Their relations with each and two alternative conceptualizations of reliability. Psychometrika. 70, 123-133.  \url{https://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf}

 
 Zinbarg, R., Yovel, I., Revelle, W. & McDonald, R. (2006).  Estimating generalizability to a universe of indicators that all have one attribute in common:  A comparison of estimators for omega.  Applied Psychological Measurement, 30, 121-144. DOI: 10.1177/0146621605278814 
 
 
}
\author{ 
 \url{https://personality-project.org/revelle.html} \cr
Maintainer: William Revelle  \email{  revelle@northwestern.edu           } 
}
\note{ omega.graph requires rgraphviz.   -- omega requires GPArotation
}
\seealso{ \code{\link{omega}},  \code{\link{make.hierarchical}}, \code{\link{ICLUST.rgraph}} }
\examples{
#24 mental tests from Holzinger-Swineford-Harman
if(require(GPArotation) ) {om24 <- omega(Harman74.cor$cov,4) } #run omega

#
#example hierarchical structure from Jensen and Weng
if(require(GPArotation) ) {jen.omega <- omega(make.hierarchical())}


}

\keyword{ multivariate }
