\name{epi}
\alias{epi}
\alias{epi.dictionary}
\docType{data}
\title{Eysenck Personality Inventory (EPI) data for 3570 participants}
\description{The EPI is and has been a very frequently administered personality test with 57 measuring two broad dimensions, Extraversion-Introversion and Stability-Neuroticism, with an additional Lie scale. Developed by Eysenck and Eysenck, 1964.  Eventually replaced with the EPQ which measures three broad dimensions.  This data set represents 3570 observations collected in the early 1990s at the Personality, Motivation and Cognition lab at Northwestern.  The data are included here as demonstration of scale construction.
}
\usage{data(epi)
data(epi.dictionary)}
\format{
  A data frame with 3570 observations on the following 57 variables.
  \describe{
    \item{\code{V1}}{a numeric vector}
    \item{\code{V2}}{a numeric vector}
    \item{\code{V3}}{a numeric vector}
    \item{\code{V4}}{a numeric vector}
    \item{\code{V5}}{a numeric vector}
    \item{\code{V6}}{a numeric vector}
    \item{\code{V7}}{a numeric vector}
    \item{\code{V8}}{a numeric vector}
    \item{\code{V9}}{a numeric vector}
    \item{\code{V10}}{a numeric vector}
    \item{\code{V11}}{a numeric vector}
    \item{\code{V12}}{a numeric vector}
    \item{\code{V13}}{a numeric vector}
    \item{\code{V14}}{a numeric vector}
    \item{\code{V15}}{a numeric vector}
    \item{\code{V16}}{a numeric vector}
    \item{\code{V17}}{a numeric vector}
    \item{\code{V18}}{a numeric vector}
    \item{\code{V19}}{a numeric vector}
    \item{\code{V20}}{a numeric vector}
    \item{\code{V21}}{a numeric vector}
    \item{\code{V22}}{a numeric vector}
    \item{\code{V23}}{a numeric vector}
    \item{\code{V24}}{a numeric vector}
    \item{\code{V25}}{a numeric vector}
    \item{\code{V26}}{a numeric vector}
    \item{\code{V27}}{a numeric vector}
    \item{\code{V28}}{a numeric vector}
    \item{\code{V29}}{a numeric vector}
    \item{\code{V30}}{a numeric vector}
    \item{\code{V31}}{a numeric vector}
    \item{\code{V32}}{a numeric vector}
    \item{\code{V33}}{a numeric vector}
    \item{\code{V34}}{a numeric vector}
    \item{\code{V35}}{a numeric vector}
    \item{\code{V36}}{a numeric vector}
    \item{\code{V37}}{a numeric vector}
    \item{\code{V38}}{a numeric vector}
    \item{\code{V39}}{a numeric vector}
    \item{\code{V40}}{a numeric vector}
    \item{\code{V41}}{a numeric vector}
    \item{\code{V42}}{a numeric vector}
    \item{\code{V43}}{a numeric vector}
    \item{\code{V44}}{a numeric vector}
    \item{\code{V45}}{a numeric vector}
    \item{\code{V46}}{a numeric vector}
    \item{\code{V47}}{a numeric vector}
    \item{\code{V48}}{a numeric vector}
    \item{\code{V49}}{a numeric vector}
    \item{\code{V50}}{a numeric vector}
    \item{\code{V51}}{a numeric vector}
    \item{\code{V52}}{a numeric vector}
    \item{\code{V53}}{a numeric vector}
    \item{\code{V54}}{a numeric vector}
    \item{\code{V55}}{a numeric vector}
    \item{\code{V56}}{a numeric vector}
    \item{\code{V57}}{a numeric vector}
  }
}
\details{
The original data were collected in a group testing framework for screening participants for subsequent studies.  The participants were enrolled in an introductory psychology class between Fall, 1991 and Spring, 1995. 

The structure of the E scale has been shown by  Rocklin and Revelle (1981) to have two subcomponents, Impulsivity and Sociability. These were subsequently used by Revelle, Humphreys, Simon and Gilliland to examine the relationship between personality, caffeine induced arousal, and cognitive performance. 
}
\source{Data from the PMC laboratory at Northwestern.  
}
\references{
Eysenck, H.J. and Eysenck, S. B.G. (1968). Manual for the Eysenck Personality Inventory.Educational and Industrial Testing Service, San Diego, CA.

Rocklin, T. and Revelle, W. (1981). The measurement of extraversion: A comparison of the Eysenck Personality Inventory and the Eysenck Personality Questionnaire.  British Journal of Social Psychology, 20(4):279-284.

}
\examples{
data(epi)
epi.keys <- make.keys(epi,list(E = c(1, 3, -5, 8, 10, 13, -15, 17, -20, 22, 25, 27,
                -29, -32, -34, -37, 39, -41, 44, 46, 49, -51, 53, 56),
   N=c(2, 4, 7, 9, 11, 14, 16, 19, 21, 23, 26, 28, 31, 33, 35, 38, 40,
     43, 45, 47, 50, 52, 55, 57),
   L = c(6, -12, -18, 24, -30, 36, -42, -48, -54),
   I =c(1, 3, -5, 8, 10, 13, 22, 39, -41), 
   S = c(-11, -15, 17, -20, 25, 27, -29, -32, -37, 44, 46, -51, 53)))
scores <- scoreItems(epi.keys,epi)
  N <- epi[abs(epi.keys[,"N"]) >0]
  E <- epi[abs(epi.keys[,"E"]) >0]
  fa.lookup(epi.keys[,1:3],epi.dictionary) #show the items and keying information
}
\keyword{datasets}
