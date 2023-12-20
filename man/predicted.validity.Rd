\name{predicted.validity}
\alias{predicted.validity}
\alias{item.validity}
\alias{validityItem}

\title{Find the predicted validities of a set of scales based on item statistics}

\description{The validity of a scale varies as a function of the number of items in the scale, their average intercorrelation, and their average validity. The asymptotic limit of a scales validity for any particular criterion is just the average validity divided by the square root of the average within scale item correlation.  \code{\link{predicted.validity}} will find the predicted validity for a set of scales (defined by a keys.list) and the average item validity for various criteria.  

The function will find (and report) scale reliabilities (using \code{\link{reliability}}) and average item validities (using \code{\link{item.validity}})
}
\usage{
predicted.validity(x, criteria, keys, scale.rel = NULL, item.val = NULL)
item.validity(x,criteria,keys) 
validityItem(x,criteria,keys)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{criteria}{Variables to predict from the scales
}
  \item{keys}{A  keys.list that defines the scales
}
  \item{scale.rel}{If not specified, these will be found.  Otherwise, this is the output from \code{\link{reliability}}.
}
  \item{item.val}{If not specified, the average item validities for each scale will be found. Otherwise use the output from \code{\link{item.validity}}
}
}
\details{
When predicting criteria from a set of items formed into scales, the validity of the scale (that is, the correlations of the scale with each criteria) is a function of the average item validity (r_y), the average intercorrelation of the items in the scale (r_x), and the number of items in the scale (n).  The limit of validity is r_y/sqrt(r_x).  

Criteria will differ in their predictability from a set of scales. These asymptotic values may be used to help the decision on which scales to develop further.  
}
\value{
\item{predicted}{The predicted validities given the scales specified}
 \item{item.validities }{The average item validities for each scale with each criterion}
\item{scale.reliabilities}{The various statistics reported by the \code{\link{reliability}} function}
 \item{asymptotic}{A matrix of the asymptotic validities}
 
}
\references{
Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 

Revelle, W. and Condon, D.M. (2019) Reliability from alpha to omega: A tutorial.  Psychological Assessment, 31, 12, 1395-1411. https://doi.org/10.1037/pas0000754.  \url{https://osf.io/preprints/psyarxiv/2y3w9} Preprint available from PsyArxiv 
}
\author{William Revelle}


\seealso{ \code{\link{reliability}}, \code{\link{scoreItems}}, \code{\link{scoreFast}}
}
\examples{
pred.bfi <- predicted.validity(psychTools::bfi[,1:25], psychTools::bfi[,26:28],
     psychTools::bfi.keys)
pred.bfi
}
\keyword{ multivariate }
\keyword{ models }