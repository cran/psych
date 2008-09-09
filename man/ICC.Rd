\name{ICC}
\alias{ICC}

\title{ Intraclass Correlations (ICC1, ICC2, ICC3 from Shrout and Fleiss) }
\description{The Intraclass correlation is used as a measure of association when studying the reliability of raters.  Shrout and Fleiss (1979) outline 6 different estimates, that depend upon the particular experimental design. All are implemented. 
}
\usage{
ICC(x,digits=2)
}

\arguments{
  \item{x}{a matrix or dataframe of ratings}
  \item{digits}{Round the output to digits}
}
\details{Shrout and Fleiss (1979) consider six cases of reliability of ratings done by k raters on n targets.

1) Each target is rated by a different set of k judges.
2) A random sample of k judges rate each target
3) A fixed set of k judges rate each target.
Then, for each of these cases, is reliability to be estimated for a single rating or for the average of k ratings?  (The 1 rating case is equivalent to the average intercorrelation, the k rating case to the Spearman Brown adjusted reliability.

}
\value{
  \item{ICC1}{Intraclass correlation for different judges}
  \item{ICC2}{Intraclass corelation for the same random judges}
  \item{ICC3}{Intraclass corelation for the same fixed judges}
  \item{ICC12}{Intraclass corelation for the average of different judges}
  \item{ICC22}{Intraclass corelation for the average of random judges}
  \item{ICC3}{Intraclass corelation for the average of fixed judges}
}
\references{ Shrout, Patrick E. and Fleiss, Joseph L. Intraclass correlations: uses in assessing rater reliability.  Psychological Bulletin, 1979, 86, 420-3428.

Revelle, W. (in prep) An introduction to psychometric theory with applications in R. Springer.  }

\author{William Revelle }

\examples{
ICC(attitude)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

