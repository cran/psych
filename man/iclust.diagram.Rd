\name{iclust.diagram}
\Rdversion{1.1}
\alias{iclust.diagram}
\alias{ICLUST.diagram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw an ICLUST hierarchical cluster structure diagram }
\description{Given a cluster structure determined by \code{\link{ICLUST}}, create a graphic structural diagram using graphic functions in the psych package To create dot code to describe the \code{\link{ICLUST}} output with more precision, use \code{\link{ICLUST.graph}}. If Rgraphviz has been successfully installed, the alternative is to use \code{\link{ICLUST.rgraph}}.
}
\usage{
iclust.diagram(ic, labels = NULL, short = FALSE, digits = 2, cex = NULL, min.size = NULL, e.size =1,colors=c("black","blue"), main = "ICLUST diagram")
}

\arguments{
  \item{ic}{Output from ICLUST}
  \item{labels}{labels for variables (if not specified as rownames in the ICLUST output}
  \item{short}{if short=TRUE, variable names are replaced with Vn}
  \item{digits}{Round the path coefficients to digits accuracy}
  \item{cex}{The standard graphic control parameter for font size modifications.  This can be used to make the labels bigger or smaller than the default values.}
  \item{min.size}{Don't provide statistics for clusters less than min.size}
  \item{e.size}{size of the ellipses with the cluster statistics.}
  \item{colors}{postive and negative }
  \item{main}{The main graphic title}
 
}
\details{iclust.diagram provides most of the power of \code{\link{ICLUST.rgraph}} without the difficulties involved in installing Rgraphviz.  It is called automatically from ICLUST. 
}
\value{
Graphical output summarizing the hierarchical cluster structure.  The graph is drawn using the diagram functions (e.g., \code{\link{dia.curve}},  \code{\link{dia.arrow}}, \code{\link{dia.rect}}, \code{\link{dia.ellipse}} ) created as a work around to Rgraphviz.
}
\references{Revelle, W. Hierarchical Cluster Analysis and the Internal Structure of Tests. Multivariate Behavioral Research, 1979, 14, 57-74.
}
\author{William Revelle
}
\note{Suggestions for improving the graphic output are welcome.   
}



\seealso{
 \code{\link{ICLUST}}
}
\examples{
v9 <- sim.hierarchical()
v9c <- ICLUST(v9)
test.data <- Harman74.cor$cov
ic.out <- ICLUST(test.data)
}

\keyword{ multivariate}
\keyword{ cluster}% __ONLY ONE__ keyword per line
\keyword{hplot}

