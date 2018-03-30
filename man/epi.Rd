\name{epi}
\alias{epi}
\alias{epi.dictionary}
\alias{epiR}
\docType{data}
\title{Eysenck Personality Inventory (EPI) data for 3570 participants}
\description{The EPI is and has been a very frequently administered personality test with 57 measuring two broad dimensions, Extraversion-Introversion and Stability-Neuroticism, with an additional Lie scale. Developed by Eysenck and Eysenck, 1964. Eventually replaced with the EPQ which measures three broad dimensions. This data set represents 3570 observations collected in the early 1990s at the Personality, Motivation and Cognition lab at Northwestern. An additional data set (epiR) has test and retest information for 474 participants. The data are included here as demonstration of scale construction and test-retest reliability.
}
\usage{data(epi)
data(epi.dictionary)
data(epiR)}
\format{
 A data frame with 3570 observations on the following 57 variables. 
 \describe{
  \item{\code{id}}{The identification number within the study}
  \item{\code{time}}{First (group testing) or 2nd time (before a lab experiment)
       for the epiR data set.}
  \item{\code{study}}{Four lab based studies and their pretest data}
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
The original data were collected in a group testing framework for screening participants for subsequent studies. The participants were enrolled in an introductory psychology class between Fall, 1991 and Spring, 1995. 

The structure of the E scale has been shown by Rocklin and Revelle (1981) to have two subcomponents, Impulsivity and Sociability. These were subsequently used by Revelle, Humphreys, Simon and Gilliland to examine the relationship between personality, caffeine induced arousal, and cognitive performance. 

The epiR data include the original group testing data and matched data for 474 participants collected several weeks later.  This is useful for showing that internal consistency estimates (e.g. \code{\link{alpha}} or \code{\link{omega}}) can be low even though the test is stable across time.  For more demonstrations of the distinction between immediate internal consistency and delayed test-retest reliability see the \code{\link{msqR}} and \code{\link{sai}} data sets and \code{\link{testRetest}}.
}
\source{Data from the PMC laboratory at Northwestern. 
}
\references{
Eysenck, H.J. and Eysenck, S. B.G. (1968). Manual for the Eysenck Personality Inventory.Educational and Industrial Testing Service, San Diego, CA.

Rocklin, T. and Revelle, W. (1981). The measurement of extraversion: A comparison of the Eysenck Personality Inventory and the Eysenck Personality Questionnaire. British Journal of Social Psychology, 20(4):279-284.

}
\examples{
data(epi)
epi.keys <- list(E = c("V1",  "V3",  "V8",  "V10", "V13", "V17", "V22", "V25", "V27", "V39",
  "V44", "V46", "V49", "V53", "V56", "-V5", "-V15", "-V20", "-V29", "-V32", "-V34","-V37",
   "-V41", "-V51"),
N = c( "V2", "V4", "V7", "V9", "V11", "V14", "V16", "V19", "V21", "V23", "V26", "V28", 
"V31", "V33", "V35", "V38", "V40","V43", "V45", "V47", "V50", "V52","V55", "V57"),
L = c("V6",  "V24", "V36", "-V12", "-V18", "-V30", "-V42", "-V48", "-V54"),
Imp = c( "V1",  "V3",  "V8",  "V10", "V13", "V22", "V39", "-V5", "-V41"),
Soc = c( "V17", "V25", "V27", "V44", "V46", "V53", "-V11", "-V15", "-V20", 
"-V29", "-V32", "-V37", "-V51")
)
scores <- scoreItems(epi.keys,epi)

keys <- make.keys(epi,epi.keys)   #the old way of making keys is to make a matrix
fa.lookup(keys[,1:3],epi.dictionary) #show the items and keying information

#a variety of demonstrations (not run) of test retest reliability versus alpha versus omega

E <- selectFromKeys(epi.keys$E)
#omega(epi[E])  #to show the low omega but high alpha of Extraversion 
#testRetest(epiR,select=E)  #test retest of the extraversion scale (.82) is higher than 
#alpha  for either the first (.77) or second administration (.74)
#Imp <- selectFromKeys(epi.keys$Imp)
#testRetest(epiR,select=Imp)  #Similarly test retest = .68 but alpha = .48 and .50.
#Soc <- selectFromKeys(epi.keys$Soc)
#testRetest(epiR,select=Soc)  # test retest =.83, alpha = .76, .75
#N <- selectFromKeys(epi.keys$N)
#testRetest(epiR, select=N) #Test retest = .8, alpha = .81, .80
}

\keyword{datasets}
