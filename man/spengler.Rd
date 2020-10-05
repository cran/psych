\name{Spengler}
\alias{Spengler}
\alias{spengler}
\alias{Damian}
\alias{Spengler.stat}
\docType{data}
\title{Project Talent data set from Marion Spengler and Rodica Damian
}
\description{Project Talent gave 440,000 US high school students a number of personality and ability tests. Of these, the data fror 346,000 were available for followup. Subsequent followups were collected 11 and 50 years later.
Marion Spengler and her colleagues Rodica Damian, and Brent Roberts reported on the stability and change across 50 years of personality and ability. Here is the correlation matrix of 25 of their variables (Spengler) as well as a slightly different set of 19 variables (Damian).
This is a nice example of mediation and regression from a correlation matrix.  (Temporarily copied from psychTools to pass CRAN checks.  )
}
\usage{data("Damian")}
\format{
 A 25 x 25 correlation matrix of demographic, personality, and ability variables, based upon 346,660 participants.
 
  
 \describe{
    \item{\code{Race/Ethnicity}}{1 = other, 2 = white/caucasian}
    \item{\code{Sex}}{1=Male, 2=Female}
    \item{\code{Age}}{Cohort =9th grade, 10th grade, 11th grade, 12th grade}
    \item{\code{Parental}}{Parental SES based upon 9 questions of home value, family income, etc.}
    \item{\code{IQ}}{Standardized composite of Verbal, Spatial and Mathematical}
    \item{\code{Sociability etc.}}{10 scales based upon prior work by Damian and Roberts}
    \item{\code{Maturity}}{A higher order factor from the prior 10 scales}
    \item{\code{Extraversion}}{The second higher order factor}
    \item{\code{Interest}}{Self reported interest in school}
    \item{\code{Reading}}{Self report reading skills}
    \item{\code{Writing}}{Self report writing skills }
    \item{\code{Responsible}}{Self reported responsibility scale}
    \item{\code{Ed.11}}{Education level at 11 year followup}
    \item{\code{Educ.50}}{Education level at 50 year followup}
    \item{\code{OccPres.11}}{Occupational Prestige at 11 year followup}
    \item{\code{OccPres.50}}{Occupational Prestige at 50 year followup}
    \item{\code{Income.11}}{Income at 11 year followup}
     \item{\code{Income.50}}{Income at 50 year followup}
    }
    }
    
\details{
Data from Project Talent was collected in 1960 on a representative sample of American high school students.  Subsequent follow up 11 and 50 years later are reported by Spengler et al (2018) and others.


}
\source{
Marion Spengler, supplementary material to Damian et al. and Spengler et al.

}
\references{
Rodica Ioana Damian and Marion Spengler and Andreea Sutu and Brent W. Roberts, 2019, Sixteen going on sixty-six: A longitudinal study of personality stability and change across 50 years 
Journal of Personality and Social Psychology. 117, (3) 274-695.

Marian Spengler and Rodica Ioana Damian and Brent W. Roberts (2018), How you behave in school predicts life success above and beyond family background, broad traits, and cognitive ability
Journal of Personality and Social Psychology, 114 (4) 600-636

}
\examples{
data(Damian)
Spengler.stat #show the basic descriptives of the original data set
psych::lowerMat(Spengler[psych::cs(IQ,Parental,Ed.11,OccPres.50),
                        psych::cs(IQ,Parental,Ed.11,OccPres.50)])
psych::setCor(OccPres.50 ~ IQ + Parental + (Ed.11),data=Spengler)
#we reduce the number of subjects for faster replication in this example
mod <- psych::mediate(OccPres.50 ~ IQ + Parental + (Ed.11),data=Spengler,
       n.iter=50,n.obs=1000) #for speed
summary(mod)

}
\keyword{datasets}
