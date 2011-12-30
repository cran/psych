\name{affect}
\alias{affect}
\alias{maps}
\alias{flat}
\docType{data}
\title{Two data sets of affect and arousal scores as a function of personality and movie conditions
}
\description{A recurring question in the study of affect is the proper dimensionality and the relationship to various personality dimensions. Here are two data sets taken from studies of mood and arousal using movies to induce affective states.  
}
\usage{data(affect)}


\details{These are data from two studies conducted in the Personality, Motivation and Cognition Laboratory at Northwestern University.  Both studies used a similar methodology: 

Collection of pretest data using 5 scales from the Eysenck Personality Inventory and items taken from the Motivational State Questionnaire (see \code{link{msq}}. In addition, state and trait anxiety measures were given.  In the ``maps" study, the Beck Depression Inventory was given also.

Then subjects were randomly assigned to one of four movie conditions: 1: Frontline.  A documentary about the liberation of the Bergen-Belsen concentration camp.  2: Halloween.  A horror film.  3: National Geographic, a nature film about the Serengeti plain.  4: Parenthood.  A comedy.  Each film clip was shown for 9 minutes.  Following this the MSQ was given again.  

Data from the MSQ were scored for Energetic and Tense Arousal (EA and TA) as well as Positive and Negative Affect (PA and NA).  

Study flat had 170 participants, study maps had 160.

These studies are described in more detail in various publications from the PMC lab.  In particular, Revelle and Anderson, 1997 and Rafaeli and Revelle (2006). 


}
\source{Data collected at the Personality, Motivation, and Cognition Laboratory, Northwestern University. 
}
\references{

William Revelle and Kristen Joan Anderson (1997) Personality, motivation and cognitive performance: Final report to the Army Research Institute on  contract MDA 903-93-K-0008

Rafaeli, Eshkol and Revelle, William (2006), A premature consensus: Are happiness and sadness truly opposite affects? Motivation and Emotion, 30, 1, 1-12.
}	
\examples{
data(affect)
describe(flat)
pairs.panels(flat[15:18],bg=c("red","black","white","blue")[maps$Film],pch=21,main="Affect varies by movies (study 'flat')")
describe(maps)
pairs.panels(maps[14:17],bg=c("red","black","white","blue")[maps$Film],pch=21,main="Affect varies by movies (study 'maps')")

}
\keyword{datasets}
