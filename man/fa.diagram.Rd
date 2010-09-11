\name{fa.diagram}
\alias{fa.graph}
\alias{fa.diagram}
\title{ Graph factor loading matrices}
\description{Factor analysis or principal components analysis results are typically interpreted in terms of the major loadings on each factor.  These structures may be represented as a table of loadings or graphically, where all loadings with an absolute value > some cut point are represented as an edge (path).
}
\usage{
fa.diagram(fa.results,Phi=NULL,sort=TRUE,labels=NULL,cut=.3,simple=TRUE,errors=FALSE,
    digits=1,e.size=.05,rsize=.15,side=2,main="Factor Analysis",cex=NULL, ...) 
fa.graph(fa.results,out.file=NULL,labels=NULL,cut=.3,simple=TRUE,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction=c("RL","TB","LR","BT"), digits=1,main="Factor Analysis",graphviz=TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fa.results}{The output of factor analysis, principal components analysis, or ICLUST analysis.  May also be a factor loading matrix from anywhere.}
  \item{Phi}{Normally not specified (it is is found in the FA, pc, or ICLUST, solution), this may be given if the input is a loadings matrix.}
  \item{out.file}{ If it exists, a dot representation of the graph will be stored here (fa.graph)}
  \item{labels}{ Variable labels }
  \item{cut}{ Loadings with abs(loading) > cut will be shown }
  \item{simple}{Only the biggest loading per item is shown}
  \item{size}{graph size }
  \item{sort}{sort the factor loadings before showing the diagram}
  \item{errors}{include error estimates (as arrows)}
  \item{e.size}{size of ellipses}
  \item{rsize}{size of rectangles}
  \item{side}{on which side should error arrows go?}
  \item{cex}{modify font size}
  \item{node.font}{what font should be used for nodes in fa.graph  }
  \item{edge.font}{what font should be used for edges in fa.graph  }
  \item{rank.direction}{ parameter passed to Rgraphviz-- which way to draw the graph }
  \item{digits}{ Number of digits to show as an edgelable }
  \item{main}{ Graphic title }
  \item{graphviz}{Should we try to use Rgraphviz for output?}
  \item{\dots}{ other parameters }
}
\details{Path diagram representations have become standard in confirmatory factor analysis, but are not yet common in exploratory factor analysis.  Representing factor structures graphically helps some people understand the structure. 

fa.diagram does not use Rgraphviz and is the preferred function.

In fa.graph, although a nice graph is drawn for the orthogonal factor case, the oblique factor drawing is acceptable, but is better if cleaned up outside of R or done using fa.diagram.

The normal input is taken from the output of either \code{\link{fa}} or \code{\link{ICLUST}}.  It is also possible to just give a factor loading matrix as input.  In this case, supplying a Phi matrix of factor correlations is also possible.

To specify the model for a structural equation confirmatory analysis of the results, use \code{\link{structure.diagram}} instead. 
}
\value{fa.diagram: A path diagram is drawn without using Rgraphviz.  This is probably the more useful function.

fa.graph: A graph is drawn using rgraphviz.  If an output file is specified, the graph instructions are also saved in the dot language.  }
 
\author{William Revelle }
\note{ fa.graph requires Rgraphviz. Because there are occasional difficulties installing Rgraphviz from Bioconductor in that some libraries are misplaced and need to be relinked, it is probably better to use fa.diagram.
}

 
\seealso{  \code{\link{omega.graph}}, \code{\link{ICLUST.graph}}, \code{\link{structure.diagram}} to convert the factor diagram to sem modeling code.
}
\examples{

test.simple <- fa(item.sim(16),2,rotate="oblimin")
#if(require(Rgraphviz)) {fa.graph(test.simple) } 
fa.diagram(test.simple)
data(bifactor)
f3 <- fa(Thurstone,3,rotate="cluster")
fa.diagram(f3,cut=.4,digits=2)
f3l <- f3$loadings
fa.diagram(f3l,main="input from a matrix")
Phi <- f3$Phi
fa.diagram(f3l,Phi=Phi,main="Input from a matrix")
fa.diagram(ICLUST(Thurstone,2,title="Two cluster solution of Thurstone"),main="Input from ICLUST")
}
\keyword{ multivariate }
\keyword{hplot}
