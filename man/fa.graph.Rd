\name{fa.graph}
\alias{fa.graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Graph factor loading matrices}
\description{Factor analysis or principal components analysis results are typically interpreted in terms of the major loadings on each factor.  These structures may be represented as a table of loadings or graphically, where all loadings with an absolute value > some cut point are represented as an edge (path).
}
\usage{
fa.graph(fa.results, out.file = NULL, labels = NULL, cut = 0.3, simple = TRUE, size = c(8, 6), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), rank.direction = "RL", digits = 1, title = "Factor Analysis", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fa.results}{The output of factor analysis or principal components analysis }
  \item{out.file}{ If it exists, a dot representation of the graph will be stored here }
  \item{labels}{ Variable labels }
  \item{cut}{ Loadings with abs(loading) > cut will be shown }
  \item{simple}{Only the biggest loading per item is shown}
  \item{size}{  }
  \item{node.font}{  }
  \item{edge.font}{ }
  \item{rank.direction}{  }
  \item{digits}{ Number of digits to show as an edgelable }
  \item{title}{ Graphic title }
  \item{\dots}{ other parameters }
}
\details{Path diagram representations have become standard in confirmatory factor analysis, but are not yet common in exploratory factor analysis.  Representing factor structures graphically helps some people understand the structure. 
}
\value{ A graph is drawn using rgraphviz.  If an output file is specified, the graph instructions are also saved in the dot language. }
  


\author{William Revelle }
\note{ Requires Rgraphviz.

As of June 1, there is an occasionally strange result when using the simple=FALSE option in Sweave.   }

 
}
\seealso{  \code{\link{omega.graph}}, \code{\link{ICLUST.graph}}}
\examples{

test.simple <- factor.pa(item.sim(16),2)
if(require(Rgraphviz)) {fa.graph(test.simple) }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

