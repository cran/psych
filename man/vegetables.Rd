\name{vegetables}
\alias{vegetables}
\alias{veg}
\docType{data}
\title{ Paired comparison of preferences for 9 vegetables}
\description{A classic data set for demonstrating Thurstonian scaling is the preference matrix of 9 vegetables from Guilford (1954).  Used by Guiford, Nunnally, and Nunally and Bernstein, this data set allows for examples of basic scaling techniques.
}
\usage{data(vegetables)}
\format{
  A data frame with 9 choices on the following 9 vegetables. The values reflect the perecentage of times where the column entry was preferred over the row entry.
  \describe{
    \item{\code{Turn}}{Turnips}
    \item{\code{Cab}}{Cabbage}
    \item{\code{Beet}}{Beets}
    \item{\code{Asp}}{Asparagus}
    \item{\code{Car}}{Carrots}
    \item{\code{Spin}}{Spinach}
    \item{\code{S.Beans}}{String Beans}
    \item{\code{Peas}}{Peas}
    \item{\code{Corn}}{Corn}
  }
}
\details{Louis L. Thurstone was a pioneer in psychometric theory and measurement of attitudes, interests, and abilities.  Among his many contributions was a systematic analysis of the process of comparative judgment (thurstone, 1927).  He considered the case of asking subjects to successively compare pairs of objects. If the same subject does this repeatedly, or if  subjects act as random replicates of each other, their judgments can be thought of as sampled from a normal distribution of underlying (latent) scale  scores for each object, Thurstone  proposed that the comparison between the value of two objects could be represented as representing the differences of the average value for each object compared to the standard deviation of the differences between objects.  The basic model is that each item has a normal distribution of response strength and that choice represents the stronger of the two response strengths.  A justification for the normality assumption is that each decision represents the sum of many independent  inputs and thus, through the central limit theorem, is normally distributed. 

Thurstone considered five different sets of assumptions about the equality and independence of the variances for each item (thurston, 1927). Torgerson expanded this analysis slightly by considering three classes of data collection (with individuals, between individuals and mixes of within and between) crossed with three sets of assumptions (equal covariance of decision process, equal correlations and small differences in variance, equal variances).  

This vegetable data set is used by Guilford and by Nunnally to demonstrate Thurstonian scaling. 
}
\source{
 Guilford, J.P. (1954)  Psychometric Methods. McGraw-Hill, New York.
}
\references{
Nunnally, J. C. (1967). Psychometric theory., McGraw-Hill, New York.\cr

Revelle, W. An introduction to psychometric theory with applications in R. (in preparation), Springer. \url{http://personality-project.org/r/book} 
}
\seealso{  \code{\link{thurstone}}}
\examples{
data(vegetables)
thurstone(veg)
}
\keyword{datasets}
