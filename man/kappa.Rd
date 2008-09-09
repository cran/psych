\name{wkappa}
\alias{wkappa}

\title{Find Cohen's kappa and weighted kappa coefficients for correlation of two raters}
\description{Cohen's kappa (Cohen, 1960) and weighted kappa (Cohen, 1968) may be used to find the agreement of two raters when using nominal scores.  

wkappa is (probability of observed matches - probability of expected matches)/(1 - probability of expected matches).  Kappa just considers the matches on the main diagonal.  Weighted kappa considers off diagonal elements as well
}
\usage{
wkappa(x, w = NULL)
}

\arguments{
  \item{x}{Either a two by n data with categorical values from 1 to p or a p x p table.  If a data array, a table will be found.}
  \item{w}{A p x p matrix of weights.  If not specified, they are set to be 1 (on the diagonal) and .5\^(distance from diagonal) off the diagonal.}
}
\details{Some categorical judgments are made using more than two outcomes.  For example, two diagnosticians might be asked to categorize patients three ways (e.g., Personality disorder, Neurosis, Psychosis).  Just as base rates affect observed cell frequencies in a two by two table, they need to be considered in the n-way table (Cohen, 1960). 
}

\value{
  \item{kappa }{Unweighted kappa}
  \item{weighted.kappa }{If weights are provided}
}
\references{Cohen, J. (1960).A coefficient of agreement for nominal scales. Educational and Psychological Measurement, 20 37-46


Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision for scaled disagreement or partial credit. Psychological Bulletin, 70, 213-220.
 }
 
\author{William Revelle }
\note{ kappa is included in psych more for completeness than necessity.  The Kappa function in the vcd package is probably preferred. 

To avoid confusion with Kappa (from vcd) or the kappa function from base, the function is named wkappa

}


\examples{
cohen <- scan()
0.44 0.05 0.01
0.07 0.20 0.03
0.09 0.05 0.06

cohen <- matrix(cohen,ncol=3)
wkappa(cohen)

fleiss <- scan()
 0.53 0.05 0.02
 0.11 0.14 0.05
 0.01 0.06 0.03

fleiss <- matrix(fleiss,ncol=3)
weights <- scan()
 1.0000 0.0000 0.4444
 0.0000 1.0000 0.6666
 0.4444 0.6666 1.0000
 
 weights <- matrix(weights,ncol=3)
 wkappa(fleiss,weights)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

