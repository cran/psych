\name{faRotations}
\alias{faRotations}
\title{Multiple rotations of factor loadings to find local minima}


\description{A dirty little secret of factor rotation algorithms is the problem of local minima (Nguyen and Waller, 2022).  Following ideas in that article, we allow for multiple random restarts and then return the global optimal solution.  Used as part of the \code{\link{fa}} function or available as a stand alone function. 
}
\usage{
faRotations(loadings, r = NULL, rotate = "oblimin", hyper = 0.15, n.rotations = 10,...)
}

\arguments{
  \item{loadings}{Factor loadings matrix from \code{\link{fa}} or \code{\link{pca}} or any N x k loadings matrix}
  \item{r}{The correlation matrix used to find the factors.  (Used to find the factor indeterminancy of the solution)}

 \item{rotate}{"none", "varimax", "Varimax","quartimax",  "bentlerT", "equamax", "varimin", "geominT" and "bifactor" are orthogonal rotations.  "Promax", "promax", "oblimin", "simplimax", "bentlerQ,  "geominQ" and "biquartimin" and "cluster" are possible oblique transformations of the solution.  Defaults to oblimin. 
 }
  \item{hyper}{The value defining when a loading is in the ``hyperplane".}

  \item{n.rotations}{The number of random restarts to use.}
  \item{...}{additional parameters, specifically, keys may be passed if using the target rotation, or delta if using geominQ, or whether to normalize if using Varimax}
}


\details{Nguyen and Waller review the problem of local minima in factor analysis.  This is a problem for all rotation algorithms, but is more so for some.  \code{\link{faRotations}} generates n.rotations different starting values and then applies the specified rotation to the original loadings using multiple start values.  Hyperplane counts and complexity indices are reported for each starting matrix, and the one with the highest hyoerplane count and the lowest complexity is returned.

Different rotation algorithms are sensitive to local minima.  The examples show this for a particularly problematic correlation matrix offered by Keith Widaman.
}

\value{
\item{loadings}{The best rotated solution}
\item{Phi}{Factor correlations}
\item{rotation.stats}{Hyperplane count, complexity.}
\item{rot.mat}{The rotation matrix used.}
}
\references{
Nguyen, H. V., & Waller, N. G. (2022, January 6). Local Minima and Factor Rotations in Exploratory Factor Analysis. Psychological Methods. Advance online publication. doi 10.1037/met0000467
}
\author{William Revelle
}
\note{Adapted from the fungible package by Waller
}


\seealso{\code{\link{fa}}, \code{\link{pca}} 
}
\examples{
f5 <- fa(bfi[,1:25],5,rotate="none")
f5r <- faRotations(f5,n.rotations=10,rotate="oblimin")   #
faRotations(f5$loadings)  #matrix input
geo <- faRotations(f5,rotate="geominQ",n.rotation=10)
 # a popular alternative, but more sensitive to local minima
describe(f5r$rotation.stats[,1:3])
describe(geo$rotation.stats[,1:3]) 

# An interesting problem from Keith Widaman
model <- 'f1 =~.8*v1+ .6*v2 + .4 *v3 
          f2 =~ .8*v4 + .6*v5 + .4*v6
          f1 ~ .5*f2'

lavParse(model, phi=TRUE)     #show the model
widaman <- sim(model)$model  
widaman     #show the correlations

test.1 <- fa(widaman, fm = "minres", nfactors = 2, rotate = "simplimax", n.rotations=1)
test.10 <- fa(widaman, fm = "minres", nfactors = 2, rotate = "simplimax")  #note the improvement
test.ob1 <- fa(widaman, fm = "minres", nfactors = 2,n.rotations=1)   #default is oblimin 

test.ob <- fa(widaman, fm = "minres", nfactors = 2)   #default is oblimin 
#stats::promax version
test.pr <- fa(widaman, fm = "minres", nfactors = 2, rotate = "promax") 
#GPArotation version
test.pro <-  fa(widaman, fm = "minres", nfactors = 2, rotate = "Promax") 
test.pro1 <- fa(widaman, fm = "minres", nfactors = 2,rotate = "Promax",n.rotations=1)  

#pca is not factor analysis, don't confuse them

test.pc1 <- pca(widaman,nfactors=2,rotate="Promax")   #defaults to one rotation
test.pc10 <- pca(widaman,nfactors=2,rotate="Promax",n.rotations=10)
test.pcob <-  pca(widaman,nfactors=2,rotate="oblimin",n.rotations=10)

#now compare the alternative solutions
#first create a data frame
phi.mat <- data.frame(rotation=cs(Truth,Simplimax.1,Simplimax.10,oblimin.1,oblimin.10,
Promax.stats,Promax.1,Promax.10, pc.Promax.1, pc.Promax.10, pc.oblimin.10),
Phi=round(c(.5,test.1$Phi[1,2],test.10$Phi[1,2],test.ob1$Phi[1,2],test.ob$Phi[1,2],
test.pr$Phi[1,2],test.pro1$Phi[1,2],test.pro$Phi[1,2],
test.pc1$Phi[1,2], test.pc10$Phi[1,2],test.pcob$Phi[1,2]),2))

#then show it
phi.mat

#Problems with Promax
pa_promax.10 <- fa(Harman74.cor$cov, 4,fm="pa" ,rotate="Promax")
pa_promax.1 <- fa(Harman74.cor$cov, 4,fm="pa" ,rotate="Promax",n.rotations=1)
round(pa_promax.10$Phi,2)
round(pa_promax.1$Phi,2)
pa_oblimin <- fa(Harman74.cor$cov, 4,fm="pa" ,rotate="oblimin")
round(pa_oblimin$Phi,2)


 }


\keyword{ multivariate }
\keyword{ models}
