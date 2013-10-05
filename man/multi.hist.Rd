\name{multi.hist}
\alias{multi.hist}
\alias{histo.density}
\title{  Multiple histograms with density and normal fits on one page}
\description{Given a matrix or data.frame, produce histograms for each variable in a "matrix" form. Include normal fits and density distributions for each plot.

The number of rows and columns may be specified, or calculated.

May be used for single variables.
}
\usage{multi.hist(x,nrow=NULL,ncol=NULL,density=TRUE,freq=FALSE,bcol="white",
      dcol=c("black","black"),dlty=c("dashed","dotted"),
      main="Histogram, Density, and Normal Fit",...)
}
\arguments{
  \item{x}{ matrix or data.frame}
  \item{nrow}{number of rows in the plot}
  \item{ncol}{number of columns in the plot}
  \item{density}{density=TRUE, show the normal fits and density distributions}
  \item{freq}{freq=FALSE shows probability densities and density distribution, freq=TRUE shows frequencies}
  \item{bcol}{Color for the bars}
  \item{dcol}{The color(s) for the normal and the density fits. Defaults to black. }
  \item{dlty}{The line type (lty) of the normal and density fits.  (specify the optional graphic parameter lwd to change the line size)}
  \item{main}{title for each panel}
  \item{...}{additional graphic parameters (e.g., col)}
}



\author{ William Revelle }

\seealso{ \code{\link{bi.bars}} for drawing pairwise histograms}
\examples{
multi.hist(sat.act) 
multi.hist(sat.act,bcol="red")
multi.hist(sat.act,dcol="blue")  #make both lines blue
multi.hist(sat.act,dcol= c("blue","red"),dlty=c("dotted", "solid")) 
multi.hist(sat.act,freq=TRUE)   #show the frequency plot
multi.hist(sat.act,nrow=2)
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
