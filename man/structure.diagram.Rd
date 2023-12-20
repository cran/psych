\name{structure.diagram}
\alias{structure.diagram}
\alias{structure.graph}
\alias{structure.sem}
\alias{lavaan.diagram}
\alias{sem.diagram}
\alias{sem.graph}



\title{Draw a structural equation model specified by two measurement models and a structural model}
\description{Graphic presentations of structural equation models are a very useful way to conceptualize sem and confirmatory factor models. Given a measurement model on x (xmodel) and on y (ymodel) as well as a path model connecting x and y (phi), draw the graph.  If the ymodel is not specified, just draw the measurement model (xmodel + phi). If the Rx or Ry matrices are specified, show the correlations between the x variables, or y variables.

Perhaps even more usefully, the function returns a model appropriate for running directly in the \emph{sem package} written by John Fox or the \emph{lavaan} package by Yves Rosseel.    For this option to work directly, it is necessary to specfy that errrors=TRUE. 

Input can be specified as matrices or the output from  \code{\link{fa}}, factanal,  or a rotation package such as \emph{GPArotation}.

For symbolic graphs, the input matrices can be character strings or mixtures of character strings and numeric vectors.

As an option, for those without Rgraphviz installed, \code{\link{structure.sem}} will just create the sem model and skip the graph. (This functionality is now included in \code{\link{structure.diagram}}.)

structure.diagram will draw the diagram without using Rgraphviz and is probably the preferred option. structure.graph will be removed eventually.

\code{\link{lavaan.diagram}} will draw either cfa or sem results from the lavaan package. It has been tested for  cfa, sem  and mimic type output.  It takes  the output object from  \emph{lavaan} and then calls \code{\link{structure.diagram}}.

}
\usage{
structure.diagram(fx, Phi=NULL,fy=NULL,labels=NULL,cut=.3,errors=FALSE,simple=TRUE,
   regression=FALSE,lr=TRUE,Rx=NULL,Ry=NULL,digits=1,e.size=.1,
    main="Structural model", ...)
structure.graph(fx,  Phi = NULL,fy = NULL, out.file = NULL, labels = NULL, cut = 0.3, 
   errors=TRUE, simple=TRUE,regression=FALSE, size = c(8, 6), 
    node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), 
    rank.direction = c("RL", "TB", "LR", "BT"), digits = 1, 
     title = "Structural model", ...)
structure.sem(fx,  Phi = NULL, fy = NULL,out.file = NULL, labels = NULL,
     cut = 0.3, errors=TRUE, simple=TRUE,regression=FALSE)
lavaan.diagram(fit,main,e.size=.1,...) 
sem.diagram(fit,main="A SEM from the sem package",...)
sem.graph(fit,out.file=NULL,main= "A SEM from the sem package",...)
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
  \item{fit}{The output from a lavaan cfa or sem}
  
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

All three structure function return a matrix of commands suitable for using in the sem or lavaan packages.  (Specify errors=TRUE to get code that will run directly in the sem package.)

The structure.graph output can be directed to an output file for post processing using the dot graphic language but requires that Rgraphviz is installed. 

lavaan.diagram will create sem, cfa, or mimic diagrams depending upon the lavaan input.

sem.diagram and sem.graph convert the output from a simple CFA done with the sem package and draw them using structure.diagram or structure.graph.
lavaan.diagram converts the output (fit) from a simple CFA done with the lavaan package and draws them using structure.diagram.
The figure is organized to show the appropriate paths between:

The correlations between the X variables (if Rx is specified) \cr
The X variables and their latent factors  (if fx is specified) \cr
The latent X and the latent Y (if Phi is specified)  \cr
The latent Y and the observed Y (if fy is specified) \cr
The correlations between the Y variables (if Ry is specified)\cr

A confirmatory factor model would specify just fx and Phi, a structural model would include fx, Phi, and fy.  The raw correlations could be shown by just including Rx and Ry.

\code{\link{lavaan.diagram}} may be called from the \code{\link{diagram}} function which also will call  \code{\link{fa.diagram}}, \code{\link{omega.diagram}} or  \code{\link{iclust.diagram}}, depending upon the class of the fit.

Other diagram functions include \code{\link{fa.diagram}}, \code{\link{omega.diagram}}.  All of these functions use the various dia functions such as \code{\link{dia.rect}}, \code{\link{dia.ellipse}}, \code{\link{dia.arrow}}, \code{\link{dia.curve}}, \code{\link{dia.curved.arrow}}, and \code{\link{dia.shape}}.
}
\value{
 \item{sem}{(invisible) a model matrix (partially) ready for input to John Fox's sem package.  It is of class ``mod" for prettier output.  }
 \item{lavaan}{(invisible) A model specification for the lavaan package.}
 \item{dotfile}{If out.file is specified, a dot language file suitable for using in a dot graphics program such as graphviz or Omnigraffle.}
 
 A graphic structural diagram in the graphics window
}

\author{William Revelle}

\seealso{  \code{\link{fa.graph}}, \code{\link{omega.graph}}, \code{\link{sim.structural}} to create artificial data sets with particular structural properties.}
\examples{

#A set of measurement and structural models
#First set up the various matrices
fx <-  matrix(c(.9,.8,.7,rep(0,9), .6,.7,-.8,rep(0,9),.5,.6,.4),ncol=3)
fy <- matrix(c(.9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)
Phi <- matrix(c(1,.35,0,0,0,
                .35,1,.5,0,0,
                0,.5, 1,0,0,
                .7,-.6, 0, 1,0,
                .0, 0, .4,0,1 ),ncol=5,byrow=TRUE)
#now draw a number of models 
f1 <- structure.diagram(fx,main = "A measurement model for x")
f2 <- structure.diagram(fx,Phi, main = "A measurement model for x")  
f3 <- structure.diagram(fy=fy, main = "A measurement model for y") 
f4 <- structure.diagram(fx,Phi,fy,main="A structural path diagram")            
f5 <- structure.diagram(fx,Phi,fy,main="A structural path diagram",errors=TRUE)

#a mimic model
fy <- matrix(c(.9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)
fx <- matrix(c(.6,.5,0,.4),ncol=2)
mimic <- structure.diagram(fx,fy=fy,simple=FALSE,errors=TRUE, main="A mimic diagram")

fy <- matrix(c(rep(.9,8),rep(0,16),rep(.8,8)),ncol=2)
structure.diagram(fx,fy=fy, e.size=.05)

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
f3 <- structure.diagram(X2,Phi=phi2,main="an alternative representation",e.size=.4)

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

###the following example is not run but is included to show how to work with lavaan
\donttest{
library(lavaan)
mod.1 <- 'A =~ A1 + A2 + A3 + A4 + A5
         C =~ C1 + C2 + C3 + C4 + C5
         E =~ E1 +E2 + E3 + E4 +E5'
fit.1 <- sem(mod.1,psychTools::bfi[complete.cases(psychTools::bfi),],std.lv=TRUE)
lavaan.diagram(fit.1)   #a normal cfa

#compare with
f3 <- fa(psychTools::bfi[complete.cases(psychTools::bfi),1:15],3)
fa.diagram(f3)


#a sem model
mod.2 <- 'A =~ A1 + A2 + A3 + A4 + A5
         C =~ C1 + C2 + C3 + C4 + C5
         E =~ E1 +E2 + E3 + E4 +E5
         E ~ A + C '        
  
fit.2 <- sem(mod.2,psychTools::bfi[complete.cases(psychTools::bfi),],std.lv=TRUE)
lavaan.diagram(fit.2, cut=0,simple=FALSE,main="sem  model")  # A SEM Model 


#a mimic model         
mod.3 <- 'A =~ A1 + A2 + A3 + A4 + A5
         C =~ C1 + C2 + C3 + C4 + C5
         E =~ E1 +E2 + E3 + E4 +E5
         A ~ age + gender
         C ~ age + gender
         E ~ age + gender'

fit.3 <- sem(mod.3,psychTools::bfi[complete.cases(psychTools::bfi),],std.lv=TRUE)
lavaan.diagram(fit.3, cut=0,simple=FALSE,main="mimic model", e.size=.03)




}
# and finally, a regression model
X7 <- matrix(c("a","b","c","d","e","f"),nrow=6)
f5 <- structure.diagram(X7,regression=TRUE,main = "Regression model")

#and a really messy regession model
x8 <- c("b1","b2","b3")
r8 <- matrix(c(1,"r12","r13","r12",1,"r23","r13","r23",1),ncol=3)
f6<- structure.diagram(x8,Phi=r8,regression=TRUE,main="Regression model")
}
\keyword{multivariate }
\keyword{hplot }
