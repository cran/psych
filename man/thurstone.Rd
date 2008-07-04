\name{thurstone}
\alias{thurstone}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Thurstone Case V scaling}
\description{Thurstone Case V scaling allows for a scaling of objects compared to other objects. As one of the cases considered by Thurstone, Case V makes the assumption of equal variances and uncorrelated distributions. }
\usage{
thurstone(x, ranks = FALSE, digits = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A square matrix or data frame of preferences, or a rectangular data frame or matrix rank order choices. }
  \item{ranks}{TRUE if rank orders are presented}
  \item{digits}{number of digits in the goodness of fit}
 
}
\details{Louis L. Thurstone was a pioneer in psychometric theory and measurement of attitudes, interests, and abilities.  Among his many contributions was a systematic analysis of the process of comparative judgment (thurstone, 1927).  He considered the case of asking subjects to successively compare pairs of objects. If the same subject does this repeatedly, or if  subjects act as random replicates of each other, their judgments can be thought of as sampled from a normal distribution of underlying (latent) scale  scores for each object, Thurstone  proposed that the comparison between the value of two objects could be represented as representing the differences of the average value for each object compared to the standard deviation of the differences between objects.  The basic model is that each item has a normal distribution of response strength and that choice represents the stronger of the two response strengths.  A justification for the normality assumption is that each decision represents the sum of many independent  inputs and thus, through the central limit theorem, is normally distributed. 

Thurstone considered five different sets of assumptions about the equality and independence of the variances for each item (thurston, 1927). Torgerson expanded this analysis slightly by considering three classes of data collection (with individuals, between individuals and mixes of within and between) crossed with three sets of assumptions (equal covariance of decision process, equal correlations and small differences in variance, equal variances).  

The data may be either a square matrix of dataframe of preferences (as proportions with the probability of the column variable being chosen over the row variable) or a matrix or dataframe of rank orders ( 1 being prefered to 2, etc.)
}
\value{
  \item{GF }{Goodness of fit 1 = 1 - sum(squared residuals/squared original) for lower off diagonal. \cr
  
    Goodness of fit 2 =  1 -  sum(squared residuals/squared original) for full matrix.}
  \item{residual }{square matrix of residuals (of class dist)}
  \item{data}{The original choice data}
  ...
}
\references{ Thurstone, L. L. (1927) A law of comparative judgments.  Psychological Review, 34, 273-286.

Revelle, W. An introduction to psychometric theory with applications in R. (in preparation), Springer. \url{http://personality-project.org/r/book} }

\author{William Revelle}

\examples{
data(vegetables)
thurstone(veg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
