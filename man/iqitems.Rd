\name{iqitems}
\alias{iqitems}
\docType{data}
\title{14 multiple choice IQ items}
\description{14 multiple choice ability items  were included as part of the Synthetic Aperture Personality Assessment (SAPA) web based personality assessment project.  The data from 1000 subjects are included here as a demonstration set for scoring multiple choice inventories and doing basic item statistics. 
}
\usage{data(iqitems)}
\format{
   A data frame with 1000 observations on the following 14 variables.
  \describe{
    \item{\code{iq1}}{In the following number series, what number comes next?  }
    \item{\code{iq8}}{Please mark the word that does not match the other words:}
    \item{\code{iq10}}{If you rearrange the letters ATNHIDLA, you will have the name of a: }
    \item{\code{iq15}}{If Jerks are Perks and some Perks are Lerks, then some Jerks are definitely Lerks. This statement is:}
    \item{\code{iq20}}{How many total legs do two ducks and three dogs have?}
    \item{\code{iq44}}{Matrix reasoning 2}
    \item{\code{iq47}}{Matrix reasoning 5}
    \item{\code{iq2}}{In the following number series, what number comes next?   1 2 4 7 12}
    \item{\code{iq11}}{The opposite of a 'stubborn' person is a '      '  person. }
    \item{\code{iq16}}{Zach is taller than Matt and Richard is shorter than Zach. Which of the following statements would be most accurate?}
    \item{\code{iq32}}{If the day before yesterday is three days after Saturday then what day is today?}
    \item{\code{iq37}}{In the following alphanumeric series, what letter comes next? Q, S, N, P, L}
    \item{\code{iq43}}{Matrix Reasoning 1}
    \item{\code{iq49}}{Matrix Reasoning 9}
  }
}
\details{14 items were sampled from 54 items given as part of the SAPA project to develop online measures of ability.  
}
\source{\url{http://personality-project.org} 
}
\references{
}
\examples{
data(iqitems)
iq.keys <- c(4,4,3,1,4,3,2,3,1,4,1,3,4,3)
score.multiple.choice(iq.keys,iqitems)
}
\keyword{datasets}
