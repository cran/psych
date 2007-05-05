\name{polychor.matrix}
\alias{polychor.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Actually, what does this do? Convert a matrix of phi coefficients to a matrix of polycoric correlations }
\description{Given a vector of a vector of frequencies, use John Fox's polycor function to convert these to polychoric correlations.   

Not ready for public consumption.
}
\usage{
polychor.matrix(x, y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix of phi coefficients }
  \item{y}{perhaps  }
}
\details{This is a stub of a function to convert a matrix of phi coefficients with an accompanying vector of frequencies into polychoric correlations.  Clearly it is not there yet, but is included in the package as a stopgap.

Please do not use.

}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\author{ William Revelle}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x,y=NULL) {
	require(polycor)
 sizex <- dim(x)[2]
if (((is.data.frame(y))|(is.matrix(y))))  sizey<-dim(y)[2] 
   else  sizey <- dim(x)[2]
        
 result<-matrix(1,nrow=sizey,ncol=sizex)   #create the output array

 xnames<- names(x)
 colnames(result)<- names(x)

if (((is.data.frame(y))|(is.matrix(y))))  rownames(result) <- names(y)
    else  rownames(result) <- names(x)
 
    
if (!((is.data.frame(y))|(is.matrix(y)))) {     #default case returns a square matrix
   for (i in 2: sizex ) {
     for (j in 1:( i-1)) {
      result[j,i]<-polychor(table(x[,j],x[,i]) )
       result[i,j] <- result[j,i]
       }
      }
      }
     else {                   #if y is input, then return the rectangular array
        for (i in 1: sizex ) {
            for (j in 1:sizey) {
                 result[j,i]<-polychor(table(x[,i],y[,j]) )
                     }
                } }
   return (result) }
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
