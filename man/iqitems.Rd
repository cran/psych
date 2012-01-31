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
     \item{\code{iq2}}{In the following number series, what number comes next?   1 2 4 7 12}
    \item{\code{iq8}}{Please mark the word that does not match the other words:}
    \item{\code{iq10}}{If you rearrange the letters ATNHIDLA, you will have the name of a: }   
     \item{\code{iq11}}{The opposite of a 'stubborn' person is a '      '  person. }
    \item{\code{iq15}}{If Jerks are Perks and some Perks are Lerks, then some Jerks are definitely Lerks. This statement is:}
       \item{\code{iq16}}{Zach is taller than Matt and Richard is shorter than Zach. Which of the following statements would be most accurate?}
    \item{\code{iq20}}{How many total legs do two ducks and three dogs have?}
    \item{\code{iq32}}{If the day before yesterday is three days after Saturday then what day is today?}
    \item{\code{iq37}}{In the following alphanumeric series, what letter comes next? Q, S, N, P, L}
    \item{\code{iq43}}{Matrix Reasoning 1}
     \item{\code{iq44}}{Matrix reasoning 2}
    \item{\code{iq47}}{Matrix reasoning 5}
    \item{\code{iq49}}{Matrix Reasoning 9}
  }
}
\details{14 items were sampled from 54 items given as part of the SAPA project  (Revelle, Wilt and Rosenthal, 2009) to develop online measures of ability. 

This data set may be used to demonstrate item response functions, \code{\link{tetrachoric}} correlations, or \code{\link{irt.fa}}.
}
\source{ The example data set is taken from the Synthetic Aperture Personality Assessment personality and ability test at \url{http://test.personality-project.org}.
} 

\references{Revelle, William, Wilt, Joshua,  and Rosenthal, Allen (2010)  Personality and Cognition: The Personality-Cognition Link. In Gruszka, Alexandra  and Matthews, Gerald   and Szymura, Blazej (Eds.) Handbook of Individual Differences in Cognition: Attention, Memory and Executive Control, Springer. 
}
\examples{
data(iqitems)
iq.keys <- c(4,4,3,1,4,3,2,3,1,4,1,3,4,3)
score.multiple.choice(iq.keys,iqitems)
#convert them to true false 
iq.scrub <- scrub(iqitems,isvalue=0)  #first get rid of the zero responses
iq.tf <-  score.multiple.choice(iq.keys,iq.scrub,score=FALSE) #convert to wrong (0) and correct (1) for analysis
describe(iq.tf)

}
\keyword{datasets}
