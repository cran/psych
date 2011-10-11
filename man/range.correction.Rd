\name{rangeCorrection}
\alias{rangeCorrection}
\title{Correct correlations for restriction of range. (Thorndike Case 2)
}
\description{In applied settings, it is typical to find a correlation between a predictor and some criterion.  Unfortunately, if the predictor is used to choose the subjects, the range of the predictor is seriously reduced.  This restricts the observed correlation to be less than would be observed in the full range of the predictor.  A correction for this problem is well known: 

Let R the unrestricted correlaton, r the restricted correlation, S the unrestricted standard deviation, s the restricted standard deviation, then

 R = rS/(s sqrt(1-r^2 + r^2(S^2/s^2)).
}
\usage{
rangeCorrection(r, sdu, sdr)
}

\arguments{
  \item{r}{The observed correlation
}
  \item{sdu}{The unrestricted standard deviation)
}
  \item{sdr}{
The restricted standard deviation
}
}
\details{
Can be used to find correlations in a restricted sample as well as the unrestricted sample.  Not the same as the correction to reliability for restriction of range.
}
\value{The corrected correlation.

}
\references{

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{http://personality-project.org/r/book/} 

Stauffer, Joseph and Mendoza, Jorge. (2001) The proper sequence for correcting correlation coefficients for range restriction and unreliability.
Psychometrika, 66, 63-68.

}
\author{
William Revelle
}

\seealso{
cRRr in the psychometric package.

}
\examples{
rangeCorrection(.33,100.32,48.19) #example from Revelle (in prep) Chapter 4.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}% __ONLY ONE__ keyword per line
