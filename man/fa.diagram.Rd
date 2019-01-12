\name{fa.diagram}
\alias{fa.graph}
\alias{fa.rgraph}
\alias{fa.diagram}
\alias{extension.diagram}
\alias{het.diagram}
\title{ Graph factor loading matrices}
\description{Factor analysis or principal components analysis results are typically interpreted in terms of the major loadings on each factor.  These structures may be represented as a table of loadings or graphically, where all loadings with an absolute value > some cut point are represented as an edge (path). \code{\link{fa.diagram}} uses the various \code{\link{diagram}} functions to draw the diagram. \code{\link{fa.graph}} generates dot code for external plotting.  \code{\link{fa.rgraph}} uses the Rgraphviz package (if available) to draw the graph. \code{\link{het.diagram}} will draw "heterarchy" diagrams of factor/scale solutions at different levels.
}
\usage{
fa.diagram(fa.results,Phi=NULL,fe.results=NULL,sort=TRUE,labels=NULL,cut=.3,
     simple=TRUE, errors=FALSE,g=FALSE,digits=1,e.size=.05,rsize=.15,side=2,
    main,cex=NULL,marg=c(.5,.5,1,.5),adj=1,ic=FALSE, ...) 
het.diagram(r,levels,cut=.3,digits=2,both=TRUE, 
        main="Heterarchy diagram",l.cex,gap.size,...) 
extension.diagram(fa.results,Phi=NULL,fe.results=NULL,sort=TRUE,labels=NULL,cut=.3,
    e.cut=.1,simple=TRUE,e.simple=FALSE,errors=FALSE,g=FALSE,
    digits=1,e.size=.05,rsize=.15,side=2,main,cex=NULL,
    marg=c(.5,.5,1,.5),adj=1,ic=FALSE, ...)

fa.graph(fa.results,out.file=NULL,labels=NULL,cut=.3,simple=TRUE,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction=c("RL","TB","LR","BT"),
     digits=1,main="Factor Analysis", ...)
fa.rgraph(fa.results,out.file=NULL,labels=NULL,cut=.3,simple=TRUE,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction=c("RL","TB","LR","BT"),
     digits=1,main="Factor Analysis",graphviz=TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fa.results}{The output of factor analysis, principal components analysis, or ICLUST analysis.  May also be a factor loading matrix from anywhere.}
  \item{Phi}{Normally not specified (it is is found in the FA, pc, or ICLUST, solution), this may be given if the input is a loadings matrix.}
  \item{fe.results}{the results of a factor extension analysis (if any)}
  \item{out.file}{ If it exists, a dot representation of the graph will be stored here (fa.graph)}
  \item{labels}{ Variable labels }
  \item{cut}{ Loadings with abs(loading) > cut will be shown }
  \item{e.cut}{extension loadings with abs(loading) > e.cut will be shown}
  \item{simple}{Only the biggest loading per item is shown}
  \item{e.simple}{Only the biggest loading per extension item is shown }
  \item{g}{Does the factor matrix reflect a g (first) factor.  If so, then draw this to the left of the variables, with the remaining factors to the right of the variables.  It is useful to turn off the simple parameter in this case.}
  \item{ic}{If drawing a cluster analysis result, should we treat it as latent variable model (ic=FALSE) or as an observed variable model (ic=TRUE) }
  \item{r}{A correlation matrix for the het.diagram function}
  \item{levels}{A list of the elements in each level}
  \item{both}{Should arrows have double heads (in het.diagram)}
  \item{size}{graph size }
  \item{sort}{sort the factor loadings before showing the diagram}
  \item{errors}{include error estimates (as arrows)}
  \item{e.size}{size of ellipses}
  \item{rsize}{size of rectangles}
  \item{side}{on which side should error arrows go?}
  \item{cex}{modify font size}
  \item{l.cex}{modify the font size in arrows, defaults to cex}
  \item{gap.size}{The gap in the arrow for the label.  Can be adjusted to compensate for variations in cex or l.cex}
  \item{marg}{sets the margins to be wider than normal, returns them to the normal size upon exit}
  \item{adj}{how many different positions (1-3) should be used for the numeric labels. Useful if they overlap each other.}
  
  \item{node.font}{what font should be used for nodes in fa.graph  }
  \item{edge.font}{what font should be used for edges in fa.graph  }
  \item{rank.direction}{ parameter passed to Rgraphviz-- which way to draw the graph }
  \item{digits}{ Number of digits to show as an edgelable }
  \item{main}{ Graphic title, defaults to "factor analyis" or "factor analysis and extension" }
  \item{graphviz}{Should we try to use Rgraphviz for output?}
  \item{\dots}{ other parameters }
}
\details{Path diagram representations have become standard in confirmatory factor analysis, but are not yet common in exploratory factor analysis.  Representing factor structures graphically helps some people understand the structure. 

By default the arrows come from the latent variables to the observed variables.  This is, of course, the factor model.  However, if the class of the object to be drawn is 'principal', then  reverse the direction of the arrows, and the 'latent' variables are no longer latent, and are shown as boxes. For cluster models, the default is to treat them as factors, but if ic =TRUE, then we treat it as a components model.

fa.diagram does not use Rgraphviz and is the preferred function.  fa.graph generates dot code to be used by an external graphics program. It does not have all the bells and whistles of fa.diagram, but these may be done in the external editor. 

Hierarchical (bifactor) models may be drawn by specifying the g parameter as TRUE.  This allows for an graphical displays of various factor transformations with a bifactor structure (e.g., \code{\link{bifactor}} and \code{\link{biquartimin}}.  See \code{\link{omega}} for an alternative way to find these structures.  

The \code{\link{het.diagram}} function will show the case of a hetarchical structure at multiple levels.  It can also be used to show the patterns of correlations between sets of scales (e.g., EPI, NEO, BFI).  The example is for showing the relationship between 3 sets of 4 variables from the Thurstone data set. The parameters l.cex and gap.size are used to adjust the font size of the labels and the gap in the lines.  

\code{\link{extension.diagram}} will draw a \code{\link{fa.extend}} result with slightly more control than using \code{\link{fa.diagram}} or the more generic \code{\link{diagram}} function.

In fa.rgraph although a nice graph is drawn for the orthogonal factor case, the oblique factor drawing is acceptable, but is better if cleaned up outside of R or done using fa.diagram. 


The normal input is taken from the output of either \code{\link{fa}} or \code{\link{ICLUST}}. This latter case displays the ICLUST results in terms of the cluster loadings, not in terms of the cluster structure.  Actually an interesting option.

It is also possible to just give a factor loading matrix as input.  In this case, supplying a Phi matrix of factor correlations is also possible. 

It is possible, using fa.graph, to export dot code for an omega solution.  fa.graph should be applied to the schmid$sl object with labels specified as the rownames of schmid$sl.  The results will need editing to make fully compatible with dot language plotting.  

To specify the model for a structural equation confirmatory analysis of the results, use \code{\link{structure.diagram}} instead. 


}
\value{fa.diagram: A path diagram is drawn without using Rgraphviz.  This is probably the more useful function.

fa.rgraph: A graph is drawn using rgraphviz.  If an output file is specified, the graph instructions are also saved in the dot language. 

fa.graph: the graph instructions are saved in the dot language. }
 
\author{William Revelle }
\note{ fa.rgraph requires Rgraphviz. Because there are occasional difficulties installing Rgraphviz from Bioconductor in that some libraries are misplaced and need to be relinked, it is probably better to use fa.diagram or fa.graph.
}

 
\seealso{  \code{\link{omega.graph}}, \code{\link{ICLUST.graph}}, \code{\link{bassAckward.diagram}}.  \code{\link{structure.diagram}} will convert the factor diagram to sem modeling code.
}
\examples{

test.simple <- fa(item.sim(16),2,rotate="oblimin")
#if(require(Rgraphviz)) {fa.graph(test.simple) } 
fa.diagram(test.simple)
f3 <- fa(Thurstone,3,rotate="cluster")
fa.diagram(f3,cut=.4,digits=2)
f3l <- f3$loadings
fa.diagram(f3l,main="input from a matrix")
Phi <- f3$Phi
fa.diagram(f3l,Phi=Phi,main="Input from a matrix")
fa.diagram(ICLUST(Thurstone,2,title="Two cluster solution of Thurstone"),main="Input from ICLUST")
het.diagram(Thurstone,levels=list(1:4,5:8,3:7))
}
\keyword{ multivariate }
\keyword{hplot}

