\name{Promax}
\alias{Promax}
\alias{TargetQ}
\alias{target.rot}
\alias{bifactor}
\alias{biquartimin}
\alias{vgQ.bimin}
\alias{vgQ.targetQ}
\title{ Perform bifactor, promax or targeted rotations and return the inter factor angles.}
\description{The bifactor rotation implements the rotation introduced by Jennrich and Bentler (2011) by calling GPForth in the GPArotation package.  promax is an oblique rotation function introduced by Hendrickson and White (1964) and implemented in the promax function in the stats package.  Unfortunately, promax does not report the inter factor correlations.  Promax does.  TargetQ does a target rotation with elements that can be missing (NA), or numeric (e.g., 0, 1).  It uses the GPArotation package. target.rot does general target rotations to an arbitrary target matrix. The default target rotation is for an independent cluster solution.
}
\usage{
bifactor(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000)
biquartimin(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000)
TargetQ(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000,Target=NULL)
Promax(x, m = 4)
target.rot(x,keys=NULL)
vgQ.bimin(L)   #called by bifactor
vgQ.targetQ(L,Target=NULL)  #called by TargetQ
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A loadings matrix}
  \item{m}{the power to which to raise the varimax loadings (for Promax)}
  \item{keys}{An arbitrary target matrix, can be composed of  any weights, but probably -1,0, 1 weights.  If missing, the target is the independent cluster structure determined by assigning every item to it's highest loaded factor.}
  \item{L}{A loadings matrix}
  \item{Target}{A matrix of values (mainly 0s, some 1s, some NAs) to which the matrix is transformed.}  
  \item{Tmat}{An initial rotation matrix}
  \item{normalize}{parameter passed to optimization routine (GPForth in the GPArotation package}
  \item{eps}{parameter passed to optimization routine (GPForth in the GPArotation package) }
  \item{maxit}{parameter passed to optimization routine (GPForth in the GPArotation package)}
}
\details{The two most useful of these four functions is probably biquartimin which implements the oblique bifactor rotation introduced by Jennrich and Bentler (2011). The second is TargetQ which allows for missing NA values in the target. Next best is the orthogonal case, bifactor.  None of these seem to be implemented in GPArotation (yet).  

The difference between biquartimin and bifactor is just that the latter is the orthogonal case which is documented in Jennrich and Bentler (2011).  It seems as if these two functions are sensitive to the starting values and random  restarts (modifying T) might be called for.

bifactor output for the 24 cognitive variable of Holzinger matches that of Jennrich and Bentler as does output for the Chen et al. problem when fm="mle" is used and the Jennrich and Bentler solution is rescaled from covariances to correlations. 

Promax is a very direct adaptation of the stats::promax function.  The addition is that it will return the interfactor correlations as well as the loadings and rotation matrix. 

In addition, it will take output from either the factanal, \code{\link{fa}} or earlier (\code{\link{factor.pa}}, \code{\link{factor.minres}} or \code{\link{principal}})  functions and select just the loadings matrix for analysis.

TargetQ implements Michael Browne's algorithm and allows specification of NA values. The Target input is a list (see examples).  It is interesting to note how powerful specifying what a factor isn't works in defining a factor.  That is, by specifying the pattern of 0s and letting most other elements be NA, the factor structure is still clearly defined.

The target.rot function is an adaptation of a function of Michael Browne's to do rotations to arbitrary target matrices.  Suggested by Pat Shrout. 

The default for target.rot is to rotate to an independent cluster structure (every items is assigned to a group with its highest loading.) 

target.rot will not handle targets that have linear dependencies (e.g., a pure bifactor model where there is a g loading and a group factor for all variables).

}
\value{
  \item{loadings }{Oblique factor loadings}
  \item{rotmat}{The rotation matrix applied to the original loadings to produce the promax soluion or the targeted matrix}
  \item{Phi}{The interfactor correlation matrix}
}
\references{Hendrickson, A. E. and  White, P. O, 1964, British Journal of Statistical Psychology, 17, 65-70.

Jennrich, Robert and Bentler, Peter (2011) Exploratory Bi-Factor Analysis. Psychometrika,  1-13

}
\author{William Revelle }
\note{A direct adaptation of the stats:promax function following suggestions to the R-help list by Ulrich Keller and John Fox. Further  modified to do  targeted rotation similar to a function of Michael Browne.

 
}
\seealso{ \code{\link{promax}}, \code{\link{factor.pa}}, \code{\link{factor.minres}}, or \code{\link{principal}} for examples of data analysis and   \code{\link{Holzinger}} or \code{\link{Bechtoldt}} for  examples of bifactor data. \code{\link{factor.rotate}} for 'hand rotation'.  }
\examples{
jen <- sim.hierarchical()
f3 <- fa(jen,3,rotate="varimax")
f3   #not a very clean solution
Promax(f3)
target.rot(f3)
m3 <- factanal(covmat=jen,factors=3)
Promax(m3)  #example of taking the output from factanal
#compare this rotation with the solution from a targeted rotation aimed for an independent cluster solution
target.rot(m3)
#now try a bifactor solution
fa(jen,3,rotate="bifactor")
fa(jen,3,rotate="biquartimin")

# compare an oblimin with a target rotation using the Browne algorithm
 #note that we are changing the factor #order (this is for demonstration only)
 Targ <- make.keys(9,list(f1=1:3,f2=7:9,f3=4:6)) 
 Targ <- scrub(Targ,isvalue=1)  #fix the 0s, allow the NAs to be estimated
 Targ <- list(Targ)  #input must be a list
#show the target
 Targ
 fa(Thurstone,3,rotate="TargetQ",Target=Targ)  #targeted rotation
#compare with oblimin
fa(Thurstone,3)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models }% __ONLY ONE__ keyword per line
