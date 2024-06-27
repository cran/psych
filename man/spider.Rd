\name{spider}
\alias{spider}
\alias{radar}

\title{Make "radar" or "spider" plots.}
\description{
Radar plots and spider plots are just two of the many ways to show multivariate data.  \code{\link{radar}}  plots correlations as vectors ranging in length from 0 (corresponding to r=-1) to 1 (corresponding to an r=1).  The vectors are arranged radially around a circle. Spider plots connect the end points of each vector. The plots are most appropriate if the variables are organized in some meaningful manner. 
}
\usage{
spider(y,x,data,labels=NULL,rescale=FALSE,center=FALSE,connect=TRUE,overlay=FALSE,
    scale=1,ncolors=31,fill=FALSE,main=NULL,...)
    
radar(x,labels=NULL,keys=NULL,center=FALSE,connect=FALSE,scale=1,ncolors=31,fill=FALSE,
    add=FALSE,linetyp="solid", main="Radar Plot",angle=0,absolute=FALSE, 
    show=TRUE,digits=2,cut=.2,circles=TRUE, shape=FALSE, clockwise=FALSE,
     delta = NULL,label.pos=NULL,position=NULL,
     xlim=c(-1,1),ylim=c(-1, 1),...) 
}

\arguments{
  \item{y}{The y variables to plot.  Each y is plotted against all the x variables}
 \item{x}{The x variables defining each line.  Each y is plotted against all the x variables}
 \item{data}{A correlation matrix from which the x and y variables are selected}
  \item{labels}{Labels (assumed to be colnames of the data matrix) for each x variable}
  \item{rescale}{If TRUE, then rescale the data to have mean 0 and sd = 1. This is used if plotting raw data rather than correlations.}
  \item{center}{if TRUE, then lines originate at the center of the plot, otherwise they start at the mid point.}
  \item{connect}{if TRUE, a spider plot is drawn, if FALSE, just a radar plot}
  \item{scale}{can be used to magnify the plot, to make small values appear larger.}
  \item{ncolors}{if ncolors > 2, then positive correlations are plotted with shades of blue and negative correlations shades of red.  This is particularly useful if fill is TRUE.  ncolors should be an odd number, so that neutral values are coded as white.  }
  \item{fill}{if TRUE, fill the polygons with colors scaled to size of correlation}
  \item{overlay}{If TRUE, plot multiple spiders on one plot, otherwise plot them as separate plots}
  \item{add}{If TRUE, add a new spider diagram to the previous one.}
  \item{linetyp}{see lty in the par options}
  \item{main}{A label or set of labels for the plots}
  \item{keys}{If a keys list is provided, then variables are grouped by the keys, with labels drawn for the key names}
  \item{angle}{Rotate the entire figure angle/nvar to the left.  Useful for drawing circumplex structures}
  \item{absolute}{If TRUE, then just use color to show correlation size}
  \item{show}{If TRUE, show the values at the end of the radar lines if they are > cut}
  \item{digits}{round the values to digits}
  \item{cut}{Just show values > cut}
  \item{circles}{Draw circles at .25, .5 and .75}
  \item{shape}{If TRUE, do not draw circles, but rather polygons with nvar sides}
  \item{clockwise}{If TRUE, organize the variables clockwise}
  \item{delta}{How far from the ends of the lines should the values be placed (defaults to 1.05 * length of line). May be vector.}
  \item{label.pos}{How far out should the labels be placed?  (defaults to 1.05 which is just outside of the outer circle.)}
  \item{position}{A way of passing the pos parameter that includes NULL as a value.  (See pos in graphics help) }
 \item{xlim}{default values may be changed for more space for labels}
 \item{ylim}{default values by be changed for more space for labelssap}
  \item{\dots}{Additional parameters can be passed to the underlying graphics call}
}
\details{Displaying multivariate profiles may be done by a series of lines (see, e.g., matplot), by colors (see, e.g., \code{\link{corPlot}}, or by radar or spider plots. Spiders are particularly suitable for showing data thought to have circumplex structure. 

To show just one variable as a function of several others, use \code{\link{radar}}.  To make multiple plots, use \code{\link{spider}}.  An additional option when comparing just a few y values is to do overlay plots.  Alternatively, set the plotting options to do several on one page.
}
\value{Either a spider or radar plot}

\author{William Revelle}

\seealso{\code{\link{corPlot}}  }

\examples{
op <- par(mfrow=c(3,2))
spider(y=1,x=2:9,data=Thurstone,connect=FALSE) #a radar plot
spider(y=1,x=2:9,data=Thurstone) #same plot as a spider plot
 spider(y=1:3,x=4:9,data=Thurstone,overlay=TRUE)
 #make a somewhat oversized plot
spider(y=26:28,x=1:25,data=cor(bfi,use="pairwise"),fill=TRUE,scale=2) 
par(op)

#another example taken from  Lippa (2001, page 193) 
lippa.df <- 
structure(list(labels = c("Assured - Dominant", "Gregarious\nExtraverted", 
"Warm\nAgreeable", "Unassuming\nIngeneous", "Unassured - Submissive", 
"Aloof\nIntroverted", "Cold\nHearted", "Arrogant\nCalculating"
), pos = c(0.8, 0.85, 0.83, 0.8, 0.75, 0.83, 0.85, 0.85), values = c(0.41, 
-0.29, -0.53, -0.61, -0.38, 0.14, 0.59, 0.6), delta = c(1.1, 
1.2, 1.2, 1.1, 1.1, 1.5, 1.2, 1.1)), row.names = c(NA, -8L), class = "data.frame")

radar(lippa.df$values,abs=TRUE,labels=lippa.df$labels,angle=90,clockwise=TRUE,lwd=3,
 label.pos=lippa.df$pos,main="Data from Lippa (2001)",scale=.9,circles=FALSE,
  cut=0,delta=lippa.df$delta)
 segments(-1,0,1,0,lwd=.2)  # Add hairline axes
 segments(0,-1,0,1,lwd=.2)
text(0,1.05,expression(italic("Masculine Instrumentality")))
text(1.05,0,expression(italic("Feminine Communion")),srt=270)

#show how to draw a hexagon
RIASEC.df <- structure(list(labels = c("Realistic", "Investigative", "Artistic", 
"Social", "Enterprising", "Conventional"), Su = c(0.84, 0.26, 
-0.35, -0.68, 0.04, -0.33), Morris = c(1.14, 0.32, -0.19, -0.38, 
0.22, 0.23)), row.names = c(NA, -6L), class = "data.frame")

 radar(RIASEC.df$Morris,RIASEC.df$labels,clockwise=TRUE,angle=0,absolute=TRUE,circl=FALSE,scale=.7,
 position=c(1,0,0,0,0,0), lwd=4,label.pos=rep(.80,6),main="",cut=0, shape=TRUE,
  delta =c(1.1,1.25,1.25, 1.25, 1.45,1.45) )
 text(-1.04,0,expression(italic("People")),srt=90)
 text(1.04,0,expression(italic("Things")),srt=270)
 text(0,.91,expression(italic("Data")))
 text(0,-.91 ,expression(italic("Ideas")))
 segments(-1,0,1,0,lwd=.2)     #add hairline axes
 segments(0,-.86,0,.86,lwd=.2)
 text(0,1.2, "Data from Su")



}
\keyword{ multivariate } 
\keyword{ hplot }
