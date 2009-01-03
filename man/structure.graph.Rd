\name{structure.graph}
\alias{structure.graph}

\title{Draw (using Graphviz) a structural equation model specified by two measurement models and a structural model}
\description{Graphic presentations of structural equation models are a very useful way to conceptualize sem and confirmatory factor models. Given a measurement model on x (xmodel) and on y (ymodel) as well as a path model connecting x and y (phi), draw the graph.  If ymodel is not specified, just draw the measurement model (xmodel + phi).

Input can be specified as matrices or the output from factanal, factor.pa, or a rotation package such as GPArotation.
}
\usage{
structure.graph(xmodel, ymodel = NULL, Phi = NULL, out.file = NULL, labels = NULL, cut = 0.3, simple = TRUE, size = c(8, 6), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), rank.direction = c("RL", "TB", "LR", "BT"), digits = 1, title = "Structural model", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmodel}{a factor model on the x variables. }
  \item{ymodel}{a factor model on the y variables (can be empty) }
  \item{Phi}{A matrix of directed relationships.  Lower diagonal values are drawn.  If the upper diagonal values match the lower diagonal, two headed arrows are drawn.  }
  \item{out.file}{name a file to send dot language instructions. }
  \item{labels}{variable labels if not specified as colnames for the matrices}
  \item{cut}{Draw paths for values > cut }
  \item{simple}{Just draw one path per x or y variable }
  \item{size}{page size of graphic }
  \item{node.font}{ font type for graph }
  \item{edge.font}{font type for graph }
  \item{rank.direction}{ Which direction should the graph be oriented }
  \item{digits}{Number of digits to draw}
  \item{title}{ Title of graphic }
  \item{\dots}{ other options to pass to Rgraphviz }
}
\details{The output can be directed to an output file for post processing using the dot graphic language.  
}
\value{
  A graphic structural diagram
}

\author{William Revelle}

\seealso{  \code{\link{fa.graph}}, \code{\link{omega.graph}}, \code{\link{sim.structural}}}
\examples{
fx <- matrix(c(.9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)
fy <- matrix(c(.6,.5,.4),ncol=1)
Phi <- matrix(c(1,0,0,0,1,0,.7,.7,1),ncol=3,byrow=TRUE)
structure.graph(fx,fy,Phi)  
}

\keyword{multivariate }
\keyword{hplot }
