\name{multi.hist}
\alias{multi.hist}
\alias{histo.density}
\alias{histBy}
\title{  Multiple histograms with density and normal fits on one page}
\description{Given a matrix or data.frame, produce histograms for each variable in a "matrix" form. Include normal fits and density distributions for each plot.

The number of rows and columns may be specified, or calculated.
May be used for single variables.
}
\usage{multi.hist(x,nrow=NULL,ncol=NULL,density=TRUE,freq=FALSE,bcol="white",
      dcol=c("black","black"),dlty=c("dashed","dotted"),
      main=NULL,mar=c(2,1,1,1),breaks=21,global=TRUE,...)
histBy(x,var,group,data=NULL,density=TRUE,alpha=.5,breaks=21,col,xlab,
            main="Histograms by group",freq=FALSE,...) 
}
\arguments{
  \item{x}{ matrix or data.frame}
  \item{var}{The variable in x to plot in histBy}
  \item{group}{The name of the variable in x to use as the grouping variable}
  \item{data}{Needs to be specified if using formula input to histBy}
  \item{nrow}{number of rows in the plot}
  \item{ncol}{number of columns in the plot}
  \item{density}{density=TRUE, show the normal fits and density distributions}
  \item{freq}{freq=FALSE shows probability densities and density distribution, freq=TRUE shows frequencies}
  \item{bcol}{Color for the bars}
  \item{dcol}{The color(s) for the normal and the density fits. Defaults to black. }
  \item{dlty}{The line type (lty) of the normal and density fits.  (specify the optional graphic parameter lwd to change the line size)}
  \item{main}{title for each panel will be set to the column name unless specified}
  \item{mar}{Specify the lower, left, upper and right hand side margin in lines -- set to be tighter than normal default of c(5,4,4,2) + .1 }
  \item{xlab}{Label for the x variable}
  \item{breaks}{The number of breaks in histBy (see hist)}
  \item{global}{If TRUE, use the same x-axis for all plots}
  \item{alpha}{The degree of transparency of the overlapping bars in histBy}
  \item{col}{A vector of colors in histBy  (defaults to the rainbow)}
  \item{...}{additional graphic parameters (e.g., col)}
}

\details{This allows for quick summaries of multiple distributions.  Particularly useful when examining the results of multiple-split halves that come from the \code{\link{reliability}} function.  

By default, will try to make a square plot with equal number of rows and columns.  However, the number of columns and rows may be specified for a particular plot.

}

\author{ William Revelle }

\seealso{ \code{\link{bi.bars}} for drawing pairwise histograms and  \code{\link{scatterHist}} for bivariate scatter and histograms}
\examples{
multi.hist(sat.act) 
multi.hist(sat.act,bcol="red")
multi.hist(sat.act,dcol="blue")  #make both lines blue
multi.hist(sat.act,dcol= c("blue","red"),dlty=c("dotted", "solid")) 
multi.hist(sat.act,freq=TRUE)   #show the frequency plot
multi.hist(sat.act,nrow=2)
histBy(sat.act,"SATQ","gender") #input by variable names
histBy(SATQ~ gender, data=sat.act) #formula input
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
