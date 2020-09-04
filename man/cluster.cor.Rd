\name{scoreOverlap}
\alias{cluster.cor}
\alias{scoreOverlap}
\title{Find correlations of composite variables (corrected for overlap) from a larger matrix.}
\description{
 Given a  n x c cluster definition matrix of -1s, 0s, and 1s (the keys) , and a n x n correlation matrix, or an N x n data matrix, find the correlations of the composite clusters.  The keys matrix can be entered by hand, copied from the clipboard (\code{\link[psychTools]{read.clipboard}}), or taken as output from the \code{\link{factor2cluster}} or \code{\link{make.keys}} functions.  Similar functionality to \code{\link{scoreItems}} which also gives item by cluster correlations.
}
\usage{
scoreOverlap(keys, r, correct = TRUE, SMC = TRUE, av.r = TRUE, item.smc = NULL, 
     impute = TRUE,select=TRUE)
cluster.cor(keys, r.mat, correct = TRUE,SMC=TRUE,item.smc=NULL,impute=TRUE)
}

\arguments{
  \item{keys}{A list of scale/cluster keys, or a  matrix of cluster keys }
  \item{r.mat}{A correlation matrix }
  \item{r}{Either a correlation matrix or a raw data matrix}
  \item{correct}{ TRUE shows both raw and corrected for attenuation correlations}
 \item{SMC}{Should squared multiple correlations be used as communality estimates for the correlation matrix? }
 \item{item.smc}{the smcs of the items may be passed into the function for speed, or calculated if SMC=TRUE }
 \item{impute}{if TRUE, impute missing scale correlations based upon the average interitem correlation, otherwise return NA.}
  \item{av.r}{Should the average r be used in correcting for overlap? smcs otherwise.}
  \item{select}{By default, just find statistics for items included in the scoring keys. This allows for finding scores from matrices with bad items if they are not included in the set of scoring keys.}
}
\details{This are two of the functions used in the SAPA (\url{https://www.sapa-project.org/}) procedures to form synthetic correlation matrices.  Given any correlation matrix of items, it is easy to find the correlation matrix of scales made up of those items. This can also be done from the original data matrix or from the correlation matrix using \code{\link{scoreItems}} which is probably preferred unless the keys are overlapping.

In the case of overlapping keys, (items being scored on multiple scales), \code{\link{scoreOverlap}} will adjust for this overlap by replacing the overlapping covariances (which are variances when overlapping) with the corresponding best estimate of an item's ``true" variance using either the average correlation or the smc estimate for that item.  This parallels the operation done when finding alpha reliability.  This is similar to ideas suggested by Cureton (1966) and Bashaw and Anderson (1966) but uses the smc or the average interitem correlation (default).

 A typical use in the SAPA project is to form item composites by clustering or factoring (see \code{\link{fa}}, \code{\link{ICLUST}}, \code{\link{principal}}), extract the clusters from these results (\code{\link{factor2cluster}}), and then form the composite correlation matrix using \code{\link{cluster.cor}}.  The variables in this reduced matrix may then be used in multiple correlatin procedures using \code{\link{mat.regress}}.
 
 The original correlation is pre and post multiplied by the (transpose) of the keys matrix. 
 
 If some correlations are missing from the original matrix this will lead to missing values (NA) for scale intercorrelations based upon those lower level correlations. If impute=TRUE (the default), a warning is issued and the correlations are imputed based upon the average correlations of the non-missing elements of each scale. 
 
 Because the alpha estimate of reliability is based upon the correlations of the items rather than upon the covariances, this estimate of alpha is sometimes called ``standardized alpha".  If the raw items are available, it is useful to compare standardized alpha with the raw alpha found using \code{\link{scoreItems}}.  They will differ substantially only if the items differ a great deal in their variances.   
 
 \code{\link{scoreOverlap}} answers an important question when developing scales and related subscales, or when comparing alternative versions of scales.  For by removing the effect of item overlap, it gives a better estimate the relationship between the latent variables estimated by the observed sum (mean) scores.
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

\references{

Bashaw, W. and Anderson Jr, H. E. (1967). A correction for replicated error in correlation coefficients. Psychometrika, 32(4):435-441.

Cureton, E. (1966). Corrected item-test correlations. Psychometrika, 31(1):93-96.

}
\author{
Maintainer: William Revelle \email{revelle@northwestern.edu}

}
\note{ See SAPA  Revelle, W., Wilt, J.,  and Rosenthal, A. (2010)  Personality and Cognition: The Personality-Cognition Link. In Gruszka, A.  and Matthews, G. and Szymura, B. (Eds.) Handbook of Individual Differences in Cognition: Attention, Memory and Executive Control, Springer. 

The second example uses the \code{\link[psychTools]{msq}} data set of 72 measures of motivational state to examine the overlap between four lower level scales and two higher level scales.
}
\seealso{ \code{\link{factor2cluster}}, \code{\link{mat.regress}}, \code{\link{alpha}}, and most importantly, \code{\link{scoreItems}}, which will do all of what cluster.cor does for most users.  cluster.cor is an important helper function for \code{\link{iclust}}
}

\examples{
#use the msq data set that shows the structure of energetic and tense arousal
small.msq <- psychTools::msq[ c("active", "energetic", "vigorous", "wakeful", 
"wide.awake", "full.of.pep", "lively", "sleepy", "tired", "drowsy","intense", 
"jittery", "fearful", "tense", "clutched.up", "quiet", "still",    "placid",
 "calm", "at.rest") ]
small.R <- cor(small.msq,use="pairwise")
keys.list <- list(
EA = c("active", "energetic", "vigorous", "wakeful", "wide.awake", "full.of.pep",
       "lively", "-sleepy", "-tired", "-drowsy"),
TA =c("intense", "jittery", "fearful", "tense", "clutched.up", "-quiet", "-still", 
       "-placid", "-calm", "-at.rest") ,

high.EA = c("active", "energetic", "vigorous", "wakeful", "wide.awake", "full.of.pep",
       "lively"),
low.EA =c("sleepy", "tired", "drowsy"),
lowTA= c("quiet", "still", "placid", "calm", "at.rest"),
highTA = c("intense", "jittery", "fearful", "tense", "clutched.up")
   ) 
   
keys <- make.keys(small.R,keys.list)
      
adjusted.scales <- scoreOverlap(keys.list,small.R)
#compare with unadjusted
confounded.scales <- cluster.cor(keys,small.R)
summary(adjusted.scales)
#note that the EA and high and low EA and TA and high and low TA 
# scale correlations are confounded
summary(confounded.scales) 
}

\keyword{ multivariate }
\keyword{ models }
