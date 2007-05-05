\name{multi.hist}
\alias{multi.hist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Multiple histograms on one screen}
\description{Given a matrix or data.frame, produce histograms for each variable in a "matrix" form. 
}
\usage{multi.hist(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ matrix or data.frame}
}

\author{ William Revelle \cr
         Northwestern University \cr
         Evanston, Illinois \cr
 \email{  revelle@northwestern.edu           } \cr
 \url{http://personality-project.org/revelle.html}}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x) {nvar <- dim(x)[2]  #number of variables
     nsize=trunc(sqrt(nvar))+1   #size of graphic
     old.par <- par(no.readonly = TRUE) # all par settings which can be changed
     par(mfrow=c(nsize,nsize))       #set new graphic parameters
     for (i in 1:nvar) {
     name=names(x)[i]                #get the names for the variables
     hist(x[,i],main=name,xlab=name) }  #draw the histograms for each variable
     on.exit(par(old.par))   #set the graphic parameters back to the original
     }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
