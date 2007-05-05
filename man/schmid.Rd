\name{schmid}
\alias{schmid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Apply the Schmid Leiman transformation to a correlation matrix}
\description{One way to find omega is to do a factor analysis of the original data set, rotate the factors obliquely, do a Schmid Leiman transformation, and then find omega. Here is the code for Schmid Leiman.  The S-L transform takes a factor or PC solution, transforms it to an oblique solution, factors the oblique solution to find a higher order (g ) factor, and then residualizes g out of the the group factors.

}
\usage{
schmid(model, nfactors = 3, pc = "pa",...)
}

\arguments{
  \item{model}{ A correlation matrix }
  \item{nfactors}{ Number of factors to extract }
  \item{pc}{ pc="pa" for principal axes, pc="pc" for principal components, pc="mle" for maximum likelihood }
  \item{...}{Allows additional parameters to be passed to the factoring routines}
}
\details{Schmid Leiman orthogonalizations are typical in the ability domain, but are not seen as often in the non-cognitive personality domain.  S-L is one way of finding the loadings of items on the general factor for estimating omega. 

A typical example would be in the study of anxiety and depression.  A general neuroticism factor (g) accounts for much of the variance, but smaller group factors of tense anxiety, panic disorder, depression, etc. also need to be considerd.

An alternative model is to consider hierarchical cluster analysis (\code{\link{ICLUST}}). 
}
\value{
  
  \item{sl }{loadings on g + nfactors group factors, communalities, uniqueness}
  \item{orthog }{original orthogonal factor loadings}
  \item{oblique}{oblique factor loadings}
  \item{fcor }{correlations among the transformed factors}
   \item{gload }{loadings of the lower order factors on g}
  ...
}
\references{\url{http://personality-project.org/r/r.omega.html} gives an example taken from Jensen and Weng, 1994 of a S-L transformation.}
\author{ William Revelle}

\seealso{  \code{\link{omega}}, \code{\link{ICLUST}},\code{\link{VSS}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, nfactors = 3, pc = "pa",...) 
{
    if (pc=="pc") {
        fact <- principal(model, nfactors,...)
    } else {if (pc=="pa") {fact <- factor.pa(model, nfactors,...) } else {
        fact <- factanal(x, covmat = model, factors = nfactors,...)
    }}
       obminfact <- oblimin(loadings(fact))
    rownames(obminfact$loadings) <- attr(model,"dimnames")[[1]]
    fload <- obminfact$loadings
    factr <- t(obminfact$Th) \%*\% (obminfact$Th)
    gfactor <- factanal(x, covmat = factr, factors = 1)
    gload <- loadings(gfactor)
    gprimaryload <- fload \%*\% gload
    colnames(gprimaryload) <- "g factor"
    u2 <- 1 - diag(fload \%*\% t(fload))
    h2 <- 1 - u2
    uniq <- 1 - fload^2
    guniq <- 1 - gprimaryload^2
    Ig <- matrix(0, ncol = nfactors, nrow = nfactors)
    diag(Ig) <- gload
    primeload <- fload \%*\% Ig
    uniq2 <- 1 - uniq - primeload^2
    sm <- sqrt(uniq2)
    schmid <- list(sl = cbind(gprimaryload, sm, h2, u2), orthog = fload, 
        fcor = factr, gloading = gload)
}

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
