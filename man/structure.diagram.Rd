\name{structure.diagram}
\alias{structure.diagram}
\alias{structure.graph}
\alias{structure.sem}


\title{Draw a structural equation model specified by two measurement models and a structural model}
\description{Graphic presentations of structural equation models are a very useful way to conceptualize sem and confirmatory factor models. Given a measurement model on x (xmodel) and on y (ymodel) as well as a path model connecting x and y (phi), draw the graph.  If ymodel is not specified, just draw the measurement model (xmodel + phi). If the Rx or Ry matrices are specified, show the correlations between the x variables, or y variables.

Perhaps even more usefully, the function returns a model appropriate for running directly in the \emph{sem package} written by John Fox.  For this option to work directly, it is necessary to specfy that errrors=TRUE. 

Input can be specified as matrices or the output from  \code{\link{fa}},  \code{\link{factor.pa}}, factanal,  or a rotation package such as \emph{GPArotation}.

For symbolic graphs, the input matrices can be character strings or mixtures of character strings and numeric vectors.

As an option, for those without Rgraphviz installed, structure.sem will just create the sem model and skip the graph. (This functionality is now included in structure.diagram.)

structure.diagram will draw the diagram without using Rgraphviz and is probably the preferred option. structure.graph will be removed eventually.

}
\usage{
structure.diagram(fx, Phi=NULL,fy=NULL,labels=NULL,cut=.3,errors=FALSE,simple=TRUE,regression=FALSE,lr=TRUE,Rx=NULL,Ry=NULL,digits=1,e.size=.1,main="Structural model", ...)
structure.graph(fx,  Phi = NULL,fy = NULL, out.file = NULL, labels = NULL, cut = 0.3, errors=TRUE, simple=TRUE,regression=FALSE, size = c(8, 6), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), rank.direction = c("RL", "TB", "LR", "BT"), digits = 1, title = "Structural model", ...)
structure.sem(fx,  Phi = NULL, fy = NULL,out.file = NULL, labels = NULL, cut = 0.3, errors=TRUE, simple=TRUE,regression=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fx}{a factor model on the x variables. }
  \item{Phi}{A matrix of directed relationships.  Lower diagonal values are drawn.  If the upper diagonal values match the lower diagonal, two headed arrows are drawn.  For a single, directed path, just the value may be specified. }
  \item{fy}{a factor model on the y variables (can be empty) }
  \item{Rx}{The correlation matrix among the x variables}
  \item{Ry}{The correlation matrix among the y variables}
  \item{out.file}{name a file to send dot language instructions. }
  \item{labels}{variable labels if not specified as colnames for the matrices}
  \item{cut}{Draw paths for values > cut }
  \item{errors}{draw an error term for observerd variables }
  \item{simple}{Just draw one path per x or y variable }
  \item{regression}{Draw a regression diagram (observed variables cause Y)}
  \item{lr}{Direction of diagram is from left to right (lr=TRUE, default) or from bottom to top (lr=FALSE) }
  \item{e.size}{size of the ellipses in structure.diagram}
  \item{main}{main title of diagram}
  \item{size}{page size of graphic }
  \item{node.font}{ font type for graph }
  \item{edge.font}{font type for graph }
  \item{rank.direction}{ Which direction should the graph be oriented }
  \item{digits}{Number of digits to draw}
  \item{title}{ Title of graphic }
  \item{\dots}{ other options to pass to Rgraphviz }
}
\details{

The recommended function is structure.diagram which does not use Rgraphviz but which does not produce dot code either.  

All three function return a matrix of commands suitable for using in the sem package.  (Specify errors=TRUE to get code that will run directly in the sem package.)

The structure.graph output can be directed to an output file for post processing using the dot graphic language but requires that Rgraphviz is installed. 

The figure is organized to show the appropriate paths between:

The correlations between the X variables (if Rx is specified) \cr
The X variables and their latent factors  (if fx is specified) \cr
The latent X and the latent Y (if Phi is specified)  \cr
The latent Y and the observed Y (if fy is specified) \cr
The correlations between the Y variables (if Ry is specified)\cr

A confirmatory factor model would specify just fx and Phi, a structural model would include fx, Phi, and fy.  The raw correlations could be shown by just including Rx and Ry.


Other diagram functions include \code{\link{fa.diagram}}, \code{\link{omega.diagram}}.  All of these functions use the various dia functions such as \code{\link{dia.rect}}, \code{\link{dia.ellipse}}, \code{\link{dia.arrow}}, \code{\link{dia.curve}}, \code{\link{dia.curved.arrow}}, and \code{\link{dia.shape}}.
}
\value{
 \item{sem}{(invisible) a model matrix (partially) ready for input to John Fox's sem package.  It is of class ``mod" for prettier output.  }
 \item{dotfile}{If out.file is specified, a dot language file suitable for using in a dot graphics program such as graphviz or Omnigraffle.}
 
 A graphic structural diagram in the graphics window
}

\author{William Revelle}

\seealso{  \code{\link{fa.graph}}, \code{\link{omega.graph}}, \code{\link{sim.structural}} to create artificial data sets with particular structural properties.}
\examples{
fx <- matrix(c(.9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)
fy <- matrix(c(.6,.5,.4),ncol=1)
Phi <- matrix(c(1,0,0,0,1,0,.7,.7,1),ncol=3,byrow=TRUE)
f1 <- structure.diagram(fx,Phi,fy,main="A structural path diagram")

#symbolic input
X2 <- matrix(c("a",0,0,"b","e1",0,0,"e2"),ncol=4)
colnames(X2) <- c("X1","X2","E1","E2")
phi2 <- diag(1,4,4)
phi2[2,1] <- phi2[1,2] <- "r"
f2 <- structure.diagram(X2,Phi=phi2,errors=FALSE,main="A symbolic model") 

#symbolic input with error 
X2 <- matrix(c("a",0,0,"b"),ncol=2)
colnames(X2) <- c("X1","X2")
phi2 <- diag(1,2,2)
phi2[2,1] <- phi2[1,2] <- "r"
f3 <- structure.diagram(X2,Phi=phi2,main="an alternative representation")

#and yet another one
X6 <- matrix(c("a","b","c",rep(0,6),"d","e","f"),nrow=6)
colnames(X6) <- c("L1","L2")
rownames(X6) <- c("x1","x2","x3","x4","x5","x6")
Y3 <- matrix(c("u","w","z"),ncol=1)
colnames(Y3) <- "Y"
rownames(Y3) <- c("y1","y2","y3")
phi21 <- matrix(c(1,0,"r1",0,1,"r2",0,0,1),ncol=3)
colnames(phi21) <- rownames(phi21) <-  c("L1","L2","Y")
f4 <- structure.diagram(X6,phi21,Y3)



# and finally, a regression model
X7 <- matrix(c("a","b","c","d","e","f"),nrow=6)
f5 <- structure.diagram(X7,regression=TRUE)

#and a really messy regession model
x8 <- c("b1","b2","b3")
r8 <- matrix(c(1,"r12","r13","r12",1,"r23","r13","r23",1),ncol=3)
f6<- structure.diagram(x8,Phi=r8,regression=TRUE)
}
\keyword{multivariate }
\keyword{hplot }
