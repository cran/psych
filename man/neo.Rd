\name{neo}
\Rdversion{1.1}
\alias{neo}
\docType{data}
\title{NEO correlation matrix from the NEO_PI_R manual}

\description{The NEO.PI.R is a widely used personality test to assess 5 broad factors (Neuroticism, Extraversion, Openness, Agreeableness and Conscientiousness) with six facet scales for each factor.  The correlation matrix of the facets is reported in the NEO.PI.R manual for 1000 subjects. 
}
\usage{data(neo)}
\format{
  A data frame of a 30 x 30 correlation matrix with the following 30 variables.
  \describe{
    \item{N1}{Anxiety}
    \item{N2}{AngryHostility}
    \item{ N3}{Depression }
    \item{ N4}{Self-Consciousness }
    \item{ N5}{Impulsiveness }
    \item{ N6}{Vulnerability }
    \item{ E1}{Warmth }
    \item{ E2}{Gregariousness }
    \item{ E3}{Assertiveness }
    \item{ E4}{Activity }
    \item{ E5}{Excitement-Seeking }
    \item{ E6}{PositiveEmotions }
    \item{ O1}{Fantasy }
    \item{ O2}{Aesthetics }
    \item{ O3}{Feelings }
    \item{ O4}{Ideas }
    \item{ O5}{Actions }
    \item{ O6}{Values }
    \item{ A1}{Trust }
    \item{ A2}{Straightforwardness }
    \item{ A3}{Altruism }
    \item{ A4}{Compliance }
    \item{ A5}{Modesty }
    \item{ A6}{Tender-Mindedness }
    \item{ C1}{Competence }
    \item{ C2}{Order }
    \item{ C3}{Dutifulness }
    \item{ C4}{AchievementStriving }
    \item{ C5}{Self-Discipline }
    \item{ C6}{Deliberation }
  }
}
\details{The past thirty years of personality research has led to a general consensus on the identification of major dimensions of personality. Variously known as the ``Big 5" or the ``Five Factor Model", the general solution represents 5 broad domains of personal and interpersonal experience.  Neuroticism and Extraversion are thought to reflect sensitivity to negative and positive cues from the environment and the tendency to withdraw or approach.  Openness is sometimes labeled as Intellect and reflects an interest in new ideas and experiences.  Agreeableness and Conscientiousness reflect tendencies to get along with others and to want to get ahead.

The factor structure of the NEO suggests five correlated factors as well as two higher level factors.  The NEO was constructed with 6 ``facets" for each of the five broad factors.  


}
\source{Costa, Paul T. and McCrae, Robert R. (1992) (NEO PI-R) professional manual. Psychological Assessment Resources, Inc. Odessa, FL. (with permission of the author and the publisher)	
}
\references{

Digman, John M. (1990) Personality structure: Emergence of the five-factor model. Annual Review of Psychology. 41, 417-440.
	

John M. Digman (1997) Higher-order factors of the Big Five. Journal of Personality and Social Psychology, 73, 1246-1256. 

McCrae, Robert R. and Costa, Paul T., Jr. (1999) A Five-Factor theory of personality. In Pervin, Lawrence A. and John, Oliver P. (eds)  Handbook of personality: Theory and research (2nd ed.) 139-153.  Guilford Press, New York. N.Y.

Revelle, William (1995), Personality processes, Annual Review of Psychology, 46, 
	295-328.

Joshua Wilt and William Revelle (2009) Extraversion and Emotional Reactivity. In Mark Leary and Rick H. Hoyle (eds). Handbook of Individual Differences in Social Behavior. Guilford Press, New York, N.Y. 

}
\examples{
data(neo)
n5 <- fa(neo,5)
neo.keys <- make.keys(30,list(N=c(1:6),E=c(7:12),O=c(13:18),A=c(19:24),C=c(25:30)))
n5p <- target.rot(n5,neo.keys) #show a targeted rotation for simple structure
n5p

}
\keyword{datasets}
