\name{fa.multi}
\alias{fa.multi}
\alias{fa.multi.diagram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multi level (hierarchical) factor analysis
}
\description{Some factor analytic solutions produce correlated factors which may in turn be factored.  If the solution has one higher order, the omega function is most appropriate.  But, in the case of multi higher order factors, then the faMulti function will do a lower level factoring and then factor the resulting correlation matrix.  Multi level factor diagrams are also shown.
}
\usage{
fa.multi(r, nfactors = 3, nfact2 = 1, n.obs = NA, n.iter = 1, rotate = "oblimin", 
 scores = "regression", residuals = FALSE, SMC = TRUE, covar = FALSE, missing = 
 FALSE,impute = "median", min.err = 0.001, max.iter = 50, symmetric = TRUE, warnings 
 =TRUE, fm = "minres", alpha = 0.1, p = 0.05, oblique.scores = FALSE, np.obs = NULL, 
 use ="pairwise", cor = "cor", ...)

fa.multi.diagram(multi.results,sort=TRUE,labels=NULL,flabels=NULL,cut=.2,gcut=.2,
    simple=TRUE,errors=FALSE,
    digits=1,e.size=.1,rsize=.15,side=3,main=NULL,cex=NULL,color.lines=TRUE
    ,marg=c(.5,.5,1.5,.5),adj=2, ...) 

}

\arguments{The arguments match those of the fa function.
  \item{r}{
A correlation matrix or raw data matrix
}
  \item{nfactors}{The desired number of factors for the lower level
}
  \item{nfact2}{The desired number of factors for the higher level
}
\item{n.obs}{Number of observations used to find the correlation matrix if using a correlation matrix.  Used for finding the goodness of fit statistics. Must be specified if using a correlaton matrix and finding confidence intervals.}
 \item{np.obs}{The pairwise number of observations. Used if using a correlation matrix and asking for a minchi solution.}
   \item{rotate}{"none", "varimax", "quartimax",  "bentlerT", "equamax", "varimin", "geominT" and "bifactor" are orthogonal rotations.  "promax", "oblimin", "simplimax", "bentlerQ,  "geominQ" and "biquartimin" and "cluster" are possible oblique transformations of the solution.  The default is to do a oblimin transformation, although  versions prior to 2009 defaulted to varimax.}
 \item{n.iter}{Number of bootstrap interations to do in fa or fa.poly}
  \item{residuals}{Should the residual matrix be shown }
  \item{scores}{the default="regression" finds factor scores using regression.  Alternatives for  estimating factor scores include simple regression ("Thurstone"), correlaton preserving ("tenBerge") as well as "Anderson" and "Bartlett" using the appropriate algorithms (see factor.scores). Although scores="tenBerge" is probably preferred for most solutions, it will lead to problems with some improper correlation matrices.  }
  \item{SMC}{Use squared multiple correlations (SMC=TRUE) or use 1  as initial communality estimate. Try using 1 if imaginary eigen values are reported. If SMC is a vector of length the number of variables, then these values are used as starting values in the case of fm='pa'. }
  \item{covar}{if covar is TRUE, factor the covariance matrix, otherwise factor the correlation matrix}
  \item{missing}{if scores are TRUE, and missing=TRUE, then impute missing values using either the median or the mean}
  \item{impute}{"median" or "mean" values are used to replace missing values}
  \item{min.err}{Iterate until the change in communalities is less than min.err}
  
  \item{max.iter}{Maximum number of iterations for convergence }
  \item{symmetric}{symmetric=TRUE forces symmetry by just looking at the lower off diagonal values}
  \item{warnings}{warnings=TRUE => warn if number of factors is too many }
  \item{fm}{factoring method  fm="minres" will do a minimum residual (OLS), fm="wls" will do a weighted least squares (WLS) solution, fm="gls" does a generalized weighted least squares (GLS), fm="pa" will do the principal factor solution, fm="ml" will do a maximum likelihood factor analysis. fm="minchi" will minimize the sample size weighted chi square when treating pairwise correlations with different number of subjects per pair.}
  \item{alpha}{alpha level for the confidence intervals for RMSEA}
  \item{p}{if doing iterations to find confidence intervals, what probability values should be found for the confidence intervals}
  \item{oblique.scores}{When factor scores are found, should they be based on the structure matrix (default) or the pattern matrix (oblique.scores=TRUE).  } 
  \item{use}{How to treat missing data, use="pairwise" is the default".  See cor for other options.}
  \item{cor}{How to find the correlations: "cor" is Pearson", "cov" is covariance, 
  "tet" is tetrachoric, "poly" is polychoric, "mixed" uses mixed cor for a mixture of tetrachorics, polychorics, Pearsons, biserials, and polyserials, Yuleb is Yulebonett, Yuleq and YuleY are the obvious Yule coefficients as appropriate}

  \item{multi.results}{The results from fa.multi}
  \item{labels}{ variable labels }
  \item{flabels}{Labels for the factors (not counting g)}
  \item{size}{size of graphics window }
  
  \item{digits}{ Precision of labels }
  \item{cex}{control font size}
  \item{color.lines}{Use black for positive, red for negative}
  \item{marg}{The margins for the figure are set to be wider than normal by default}
  \item{adj}{Adjust the location of the factor loadings to vary as factor mod 4 + 1}
 
    \item{main}{ main figure caption }
  \item{\dots}{additional parameters, specifically, keys may be passed if using the target rotation, or delta if using geominQ, or whether to normalize if using Varimax.  In addition, for fa.multi.diagram, other options to pass into the graphics packages }
  \item{e.size}{the size to draw the ellipses for the factors. This is scaled by the  number of variables.}
  \item{cut}{Minimum path coefficient to draw}
  \item{gcut}{Minimum general factor path to draw}
  \item{simple}{draw just one path per item}
  \item{sort}{sort the solution before making the diagram}
  \item{side}{on which side should errors  be drawn?}
  \item{errors}{show the error estimates}
  \item{rsize}{size of the rectangles} 
}
\details{
See \code{\link{fa}} and \code{\link{omega}} for a discussion of factor analysis and of the case of one higher order factor.
}
\value{
\item{f1}{The standard output from a factor analysis from  \code{\link{fa}} for the raw variables}
\item{f2}{The standard output from a factor analysis from  \code{\link{fa}} for the correlation matrix of the level 1 solution. }
}

\references{
Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{http://personality-project.org/r/book/} 
}
\author{
William Revelle
}
\note{
This is clearly an early implementation (Feb 14 2016) which might be improved.
}



\seealso{
\code{\link{fa}}, \code{\link{omega}}
}
\examples{
f31 <- fa.multi(Thurstone,3,1) #compare with \code{\link{omega}}
f31
fa.multi.diagram(f31)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}
