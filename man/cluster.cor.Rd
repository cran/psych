\name{cluster.cor}
\alias{cluster.cor}

\title{Find correlations of composite variables from a larger matrix}
\description{
 Given a  n x c cluster definition matrix of -1s, 0s, and 1s (the keys) , and a n x n correlation matrix, find the correlations of the composite clusters.  The keys matrix can be entered by hand, copied from the clipboard (\code{\link{read.clipboard}}), or taken as output from the \code{\link{factor2cluster}} function. 
}
\usage{
cluster.cor(keys, r.mat, correct = TRUE,digits=2,SMC=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keys}{ A matrix of cluster keys }
  \item{r.mat}{A correlation matrix }
  \item{correct}{ TRUE shows both raw and corrected for attenuation correlations}
 \item{digits}{round off answer to digits}
 \item{SMC}{Should squared multiple correlations be used as communality estimates for the correlation matrix? }
}
\details{This is one of the functions used in the SAPA procedures to form synthetic correlation matrices.  Given any correlation matrix of items, it is easy to find the correlation matrix of scales made up of those items. This can also be done from the original data matrix using \code{\link{score.items}}.

 A typical use in the SAPA project is to form item composites by clustering or factoring (see \code{\link{factor.pa}}, \code{\link{ICLUST}}, \code{\link{principal}}), extract the clusters from these results (\code{\link{factor2cluster}}), and then form the composite correlation matrix using \code{\link{cluster.cor}}.  The variables in this reduced matrix may then be used in multiple correlatin procedures using \code{\link{mat.regress}}.
 
 The original correlation is pre and post multiplied by the (transpose) of the keys matrix. 
 
 If some correlations are missing from the original matrix this will lead to missing values (NA) for scale intercorrelations based upon those lower level correlations.  
 
 Because the alpha estimate of reliability is based upon the correlations of the items rather than upon the covariances, this estimate of alpha is sometimes called ``standardized alpha".  If the raw items are available, it is useful to compare standardized alpha with the raw alpha found using \code{\link{score.items}}.  They will differ substantially only if the items differ a great deal in their variances.   
}
\value{
  \item{cor }{the (raw) correlation matrix of the clusters}
  \item{sd }{standard deviation of the cluster scores}
  \item{corrected }{raw correlations below the diagonal, alphas on diagonal, disattenuated above diagonal}
  \item{alpha}{The (standardized) alpha reliability of each scale.}
  \item{G6}{Guttman's Lambda 6 reliability estimate is based upon the smcs for each item in a scale.  G6 uses the smc based upon the entire item domain.}
  \item{av.r}{The average inter item correlation within a scale}
 \item{size}{How many items are in each cluster?}
 
}

\author{
Maintainer: William Revelle \email{revelle@northwestern.edu}

}
\note{ See SAPA  Revelle, W., Wilt, J.,  and Rosenthal, A. (2009)  Personality and Cognition: The Personality-Cognition Link. In Gruszka, A.  and Matthews, G. and Szymura, B. (Eds.) Handbook of Individual Differences in Cognition: Attention, Memory and Executive Control, Springer. 
}
\seealso{ \code{\link{factor2cluster}}, \code{\link{mat.regress}}, \code{\link{alpha.scale}}, \code{\link{score.items}}}
\examples{
\dontrun{
data(attitude)
keys <- matrix(c(1,1,1,0,0,0,0,
                 0,0,0,1,1,1,1),ncol=2)
colnames(keys) <- c("first","second")
r.mat <- cor(attitude)
cluster.cor(keys,r.mat)
}
#$cor
#       first second
#first    1.0    0.6
#second   0.6    1.0
#
#$sd
# first second 
#  2.57   3.01 
#
#$corrected
#       first second
#first   0.82   0.77
#second  0.60   0.74
#
#$size
# first second 
#     3      4 


}

\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
