\name{schmid}
\alias{schmid}
\title{Apply the Schmid Leiman transformation to a correlation matrix}
\description{One way to find omega is to do a factor analysis of the original data set, rotate the factors obliquely, do a Schmid Leiman transformation, and then find omega. Here is the code for Schmid Leiman.  The S-L transform takes a factor or PC solution, transforms it to an oblique solution, factors the oblique solution to find a higher order (g ) factor, and then residualizes g out of the the group factors.

}
\usage{
schmid(model, nfactors = 3, fm = "minres",digits=2,rotate="oblimin",n.obs=NA,option="equal",...)
}

\arguments{
  \item{model}{ A correlation matrix }
  \item{nfactors}{ Number of factors to extract }
  \item{fm}{ fm="pa" for principal axes, fm="pc" for principal components, fm = "minres" for minimum residual (OLS), pc="mle" for maximum likelihood }
  \item{digits}{if digits not equal NULL, rounds to digits}
  \item{rotate}{The default, oblimin, produces somewhat more correlated factors than the alternative, simplimax. The third option is the promax criterion}
  \item{n.obs}{Number of observations, used to find fit statistics if specified.  Will be calculated if input is raw data}
  \item{option}{When asking for just two group factors, option can be for "equal", "first" or "second"}
  \item{...}{Allows additional parameters to be passed to the factoring routines}
}
\details{Schmid Leiman orthogonalizations are typical in the ability domain, but are not seen as often in the non-cognitive personality domain.  S-L is one way of finding the loadings of items on the general factor for estimating omega. 

A typical example would be in the study of anxiety and depression.  A general neuroticism factor (g) accounts for much of the variance, but smaller group factors of tense anxiety, panic disorder, depression, etc. also need to be considerd.

An alternative model is to consider hierarchical cluster analysis techniques such as \code{\link{ICLUST}}.

Requires the GPArotation package.

Although 3 factors are the minimum number necessary to define the solution uniquely, it is occasionally useful to allow for a two factor solution.  There are three possible options for this condition: setting the general factor loadings between the two lower order factors to be "equal" which will be the sqrt(oblique correlations between the factors) or to "first" or "second" in which case the general factor is equated with either the first or second group factor. A  message is issued suggesting that the model is not really well defined. 



}
\value{
  
  \item{sl }{loadings on g + nfactors group factors, communalities, uniqueness}
  \item{orthog }{original orthogonal factor loadings}
  \item{oblique}{oblique factor loadings}
  \item{phi }{correlations among the transformed factors}
   \item{gload }{loadings of the lower order factors on g}
  ...
}
\references{\url{http://personality-project.org/r/r.omega.html} gives an example taken from Jensen and Weng, 1994 of a S-L transformation.}
\author{ William Revelle}

\seealso{  \code{\link{omega}}, \code{\link{omega.graph}}, \code{\link{fa.graph}}, \code{\link{ICLUST}},\code{\link{VSS}}}
\examples{
jen <- sim.hierarchical()  #create a hierarchical demo
if(!require(GPA.rotation)) {message("I am sorry, you must have GPArotation installed to use schmid.")} else {jen <- schmid(jen,digits=2)   #use the oblimin rotation
p.jen <- schmid(jen,rotate="promax") #use the promax rotation
}
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
