\name{iclust.diagram}
\Rdversion{1.1}
\alias{iclust.diagram}
\alias{ICLUST.diagram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw an ICLUST hierarchical cluster structure diagram }
\description{Given a cluster structure determined by \code{\link{ICLUST}}, create a graphic structural diagram using graphic functions in the psych package To create dot code to describe the \code{\link{ICLUST}} output with more precision, use \code{\link{ICLUST.graph}}. If Rgraphviz has been successfully installed, the alternative is to use \code{\link{ICLUST.rgraph}}.
}
\usage{
iclust.diagram(ic, labels = NULL, short = FALSE, digits = 2, cex = NULL, 
    min.size = NULL,
     e.size =1,
     colors=c("black","blue"), 
     main = "ICLUST diagram",
     cluster.names=NULL,marg=c(.5,.5,1.5,.5),plot=TRUE, bottomup=TRUE)
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
  \item{cluster.names}{Normally, clusters are named sequentially C1 ... Cn.  If cluster.names are specified, then these values will be used instead.}
  \item{marg}{Sets the margins to be narrower than the default values. Resets them upon return}
  \item{plot}{If plot is TRUE, then draw the diagram,  if FALSE, then just return the veriable order from the plot}
  \item{bottomup}{Which way to draw the arrows.  TRUE means from the items to the clusters.  See note. }
 
}
\details{iclust.diagram provides most of the power of \code{\link{ICLUST.rgraph}} without the difficulties involved in installing Rgraphviz.  It is called automatically from ICLUST. 

Following a request by Michael Kubovy, cluster.names may be specified to replace the normal C1 ... Cn names.  

If access to a dot language graphics program is available, it is probably better to use the iclust.graph function to get dot output for offline editing.

Until 3/11/2 arrows went from clusters to items.  The default value for bottomup has been changed to draw from items to clusters.  To draw the old way, set bottomup=TRUE.}
\value{
Graphical output summarizing the hierarchical cluster structure.  The graph is drawn using the diagram functions (e.g., \code{\link{dia.curve}},  \code{\link{dia.arrow}}, \code{\link{dia.rect}}, \code{\link{dia.ellipse}} ) created as a work around to Rgraphviz.

Also returned (invisibly) is a vector of variable names ordered by their location in the tree diagram. The plot option suppresses the plot for speed. 
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
#now show how to relabel clusters
ic.bfi <- iclust(psychTools::bfi[1:25],beta=3) #find the clusters
cluster.names <- rownames(ic.bfi$results) #get the old names
#change the names to the desired ones
cluster.names[c(16,19,18,15,20)] <- c("Neuroticism","Extra-Open","Agreeableness",
      "Conscientiousness","Open")
#now show the new names
iclust.diagram(ic.bfi,cluster.names=cluster.names,min.size=4,e.size=1.75)

}

\keyword{ multivariate}
\keyword{ cluster}% __ONLY ONE__ keyword per line
\keyword{hplot}

