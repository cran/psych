\name{fa.graph}
\alias{fa.graph}
\alias{fa.diagram}
\title{ Graph factor loading matrices}
\description{Factor analysis or principal components analysis results are typically interpreted in terms of the major loadings on each factor.  These structures may be represented as a table of loadings or graphically, where all loadings with an absolute value > some cut point are represented as an edge (path).
}
\usage{
fa.diagram(fa.results,sort=TRUE,labels=NULL,cut=.3,simple=TRUE,errors=FALSE,
    digits=1,e.size=.05,rsize=.15,side=2,main="Factor Analysis", ...) 
fa.graph(fa.results,out.file=NULL,labels=NULL,cut=.3,simple=TRUE,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction=c("RL","TB","LR","BT"), digits=1,main="Factor Analysis",graphviz=TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fa.results}{The output of factor analysis or principal components analysis }
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
  \item{node.font}{  }
  \item{edge.font}{ }
  \item{rank.direction}{  }
  \item{digits}{ Number of digits to show as an edgelable }
  \item{main}{ Graphic title }
  \item{graphviz}{Should we try to use Rgraphviz for output?}
  \item{\dots}{ other parameters }
}
\details{Path diagram representations have become standard in confirmatory factor analysis, but are not yet common in exploratory factor analysis.  Representing factor structures graphically helps some people understand the structure. 

fa.diagram does not use Rgraphviz and is the preferred function.


In fa.graph, although a nice graph is drawn for the orthogonal factor case, the oblique factor drawing is acceptable, but is better if cleaned up outside of R.
}
\value{fa.diagram: A path diagram is drawn without using Rgraphviz.  This is probably the more useful function.

fa.graph: A graph is drawn using rgraphviz.  If an output file is specified, the graph instructions are also saved in the dot language.  }
  


\author{William Revelle }
\note{ fa.graph requires Rgraphviz. Because there are occasional difficulties installing Rgraphviz from Bioconductor in that some libraries are misplaced and need to be relinked, it is probably better to use fa.diagram.

As of June 1, 2007 there is an occasionally strange result when using the simple=FALSE option in Sweave.   }

 
\seealso{  \code{\link{omega.graph}}, \code{\link{ICLUST.graph}}}
\examples{


test.simple <- fa(item.sim(16),2,rotate="oblimin")
if(require(Rgraphviz)) {fa.graph(test.simple) }
fa.diagram(test.simple)
data(bifactor)
fa.diagram(fa(Thurstone,3,rotate="cluster"),cut=.4,digits=2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{hplot}

