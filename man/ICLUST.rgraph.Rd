\name{ICLUST.rgraph}
\alias{ICLUST.rgraph}

\title{ Draw an ICLUST graph using the Rgraphviz package }
\description{Given a cluster structure determined by \code{\link{ICLUST}}, create a rgraphic directly using Rgraphviz.  To create dot code to describe the \code{\link{ICLUST}} output with more precision, use \code{\link{ICLUST.graph}}.  As an option, dot code is also generated and saved in a file. To  use the dot code, use either \url{http://www.graphviz.org/} Graphviz or a commercial viewer (e.g., OmniGraffle).
}
\usage{
ICLUST.rgraph(ic.results, out.file = NULL, min.size = 1, short = FALSE, labels = NULL, size = c(8, 6), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), rank.direction=c("RL","TB","LR","BT"), digits = 2, title = "ICLUST",label.font=2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
   \item{ic.results}{output list from ICLUST }
   \item{out.file}{ File name to save optional dot code. }
  \item{min.size}{draw a smaller node (without all the information) for clusters < min.size -- useful for large problems}
  \item{short}{if short==TRUE, don't use variable names}
  \item{labels}{vector of text labels (contents) for the variables}
  \item{size}{size of output }
  \item{node.font}{ Font to use for nodes in the graph }
  \item{edge.font}{ Font to use for the labels of the arrows (edges)}
  \item{rank.direction}{LR or TB  or RL  }
  \item{digits}{ number of digits to show }
  \item{title}{ any title }
  \item{label.font}{The variable labels can be a different size than the other nodes.  This is particularly helpful if the number of variables is large or the labels are long.}
  \item{\dots}{ other options to pass }
}
\details{
  Will create (or overwrite) an output file and print out the dot code to show a cluster structure. This dot file may be imported directly into a dot viewer (e.g.,  \url{http://www.graphviz.org/}).  The "dot" language is a powerful graphic description language that is particulary appropriate for viewing cluster output.  Commercial graphics programs (e.g., OmniGraffle) can also read (and clean up) dot files.  
  
ICLUST.rgraph takes the output from \code{\link{ICLUST}} results and processes it to provide a pretty picture of the results.  Original variables shown as rectangles and ordered on the left hand side (if rank direction is RL) of the graph.  Clusters are drawn as ellipses and include the alpha, beta, and size of the cluster.  Edges show the cluster intercorrelations.

It is possible to trim the output to not show all cluster information. Clusters < min.size are shown as small ovals without alpha, beta, and size information.
}
\value{Output is a set of dot commands written either to console or to the output file.  These commands may then be used as input to any "dot" viewer, e.g., Graphviz.

ICLUST.rgraph is a version of \code{\link{ICLUST.graph}} that uses Rgraphviz to draw on the screen as well.  

Additional output is drawn to main graphics screen. }
\references{ ICLUST: \url{http://personality-project.org/r/r.iclust.html}}
\author{
 \email{revelle@northwestern.edu } \cr
 \url{http://personality-project.org/revelle.html}}

\seealso{  \code{\link{VSS.plot}},  \code{\link{ICLUST}}}
\note{ Requires Rgraphviz}


\examples{
test.data <- Harman74.cor$cov

if(require(Rgraphviz) ) {ic.out <- ICLUST(test.data,labels=colnames(test.data)) }  

}
\keyword{ multivariate}
\keyword{ cluster}% __ONLY ONE__ keyword per line
\keyword{hplot}
