\name{error.crosses}
\alias{error.crosses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot x and y error bars }
\description{Given two vectors of data, plot the means and show standard errors in both X and Y directions. 
}
\usage{
error.crosses(x, y, labels = NULL, pos = NULL, arrow.len = 0.2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A vector of summary statistics (from Describe) }
  \item{y}{ A second vector of summary statistics (also from Describe)}
  \item{labels}{ name the pair }
  \item{pos}{Labels are located where with respect to the mean?  }
  \item{arrow.len}{ Arrow length }
  \item{\dots}{ Other parameters for plot  }
}
\details{For an example of two way error bars describing the effects of mood manipulations upon positive and negative affect, see \url{http://personality-project.org/revelle/publications/happy-sad-appendix/FIG.A-6.pdf})
}
}
\author{ William Revelle \cr     
 \email{revelle@northwestern.edu}
}

\examples{

#desc <- describe(attitude)
#x <- desc[1,]
#y <- desc[2,]
#plot(x$mean,y$mean)   #in graphics window
#error.crosses(x,y)    #in graphics window

## The function is currently defined as

function (x,y,labels=NULL,pos=NULL,arrow.len=.2,...)  # x  and y are data frame with 
    {z <- dim(x)[1]
     if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     if (length(labels)==0) lab <- rep("",z) else lab <-labels
        for (i in 1:z)  
    	{xcen <- x$mean[i]
    	 ycen <- y$mean[i]
    	 xse  <- x$se[i]
    	 yse <-  y$se[i]
    	 arrows(xcen-xse,ycen,xcen+xse,ycen,length=arrow.len, angle = 90,
    	 code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 arrows(xcen,ycen-yse,xcen,ycen+yse,length=arrow.len, angle = 90,
    	 code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	text(xcen,ycen,labels=lab[i],pos=pos[i],cex=1,offset=arrow.len+1) 
    	#puts in labels for all points
    	}	
   }

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
