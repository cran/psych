\name{bfi}
\alias{bfi}
\docType{data}
\title{25 Personality items representing 5 factors}
\description{25 personality self report items taken from the International Personality Item Pool (ipip.ori.org) were included as part of the Synthetic Aperture Personality Assessment (SAPA) web based personality assessment project.  The data from 2800 subjects are included here as a demonstration set for scale construction, factor analysis, and Item Response Theory analysis.  Three additional demographic variables (sex, education, and age) are also included.
}
\usage{data(bfi)}
\format{
  A data frame with 2800 observations on the following 28 variables. (The q numbers are the SAPA item numbers).
  \describe{
    \item{\code{A1}}{Am indifferent to the feelings of others. (q_146)}
    \item{\code{A2}}{Inquire about others' well-being. (q_1162)}
    \item{\code{A3}}{Know how to comfort others. (q_1206) }
    \item{\code{A4}}{Love children. (q_1364)}
    \item{\code{A5}}{Make people feel at ease. (q_1419)}
    \item{\code{C1}}{Am exacting in my work. (q_124)}
    \item{\code{C2}}{Continue until everything is perfect. (q_530)}
    \item{\code{C3}}{Do things according to a plan. (q_619)}
    \item{\code{C4}}{Do things in a half-way manner. (q_626)}
    \item{\code{C5}}{Waste my time. (q_1949)}
    \item{\code{E1}}{Don't talk a lot. (q_712)}
    \item{\code{E2}}{Find it difficult to approach others. (q_901)}
    \item{\code{E3}}{Know how to captivate people. (q_1205)}
    \item{\code{E4}}{Make friends easily. (q_1410)}
    \item{\code{E5}}{Take charge. (q_1768)}
    \item{\code{N1}}{Get angry easily. (q_952)}
    \item{\code{N2}}{Get irritated easily. (q_974)}
    \item{\code{N3}}{Have frequent mood swings. (q_1099}
    \item{\code{N4}}{Often feel blue. (q_1479)}
    \item{\code{N5}}{Panic easily. (q_1505)}
    \item{\code{O1}}{Am full of ideas. (q_128)}
    \item{\code{O2}}{Avoid difficult reading material.(q_316)}
    \item{\code{O3}}{Carry the conversation to a higher level. (q_492)}
    \item{\code{O4}}{Spend time reflecting on things. (q_1738)}
    \item{\code{O5}}{Will not probe deeply into a subject. (q_1964)}
   \item{\code{gender}}{Males = 1, Females =2}
   \item{\code{education}}{1 = HS, 2 = finished HS, 3 = some college, 4 = college graduate 5 = graduate degree}
   \item{\code{age}}{age in years}
}
}
\details{The first 25 items are organized by five putative factors: Agreeableness, Conscientiousness, Extraversion, Neuroticism, and Opennness.  The scoring key is created using  \code{\link{make.keys}}, the scores are found using  \code{\link{score.items}}.

These five factors are a useful example of using \code{\link{irt.fa}} to do Item Response Theory based latent factor analysis of the \code{\link{polychoric}} correlation matrix.  The endorsement plots for each item, as well as the item information functions reveal that the items differ in their quality.

The item data were collected using a 6 point response scale: 
1 Very Inaccurate
2 Moderately Inaccurate
3 Slightly Inaccurate
4 Slightly Accurate
5 Moderately Accurate
6 Very Accurate

as part of the Synthetic Apeture Personality Assessment (SAPA) project.  To see an example of the data collection technique, visit \url{http://test.personality-project.org}.  The items given were sampled from the International Personality Item Pool of Lewis Goldberg using the sampling technique of SAPA.  This is a sample data set taken from the much larger SAPA data bank.

}
\source{The items are from the ipip (Goldberg, 1999).  The data are from the SAPA project (Revelle, Wilt and Rosenthal, 2010) , collected Spring, 2010.
}
\references{Goldberg, L.R. (1999) A broad-bandwidth, public domain, personality inventory measuring the lower-level facets of several five-factor models. In Mervielde, I. and Deary, I. and De Fruyt, F. and Ostendorf, F. (eds) Personality psychology in Europe. 7. Tilburg University Press. Tilburg, The Netherlands.

Revelle, W., Wilt, J.,  and Rosenthal, A. (2010)  Personality and Cognition: The Personality-Cognition Link. In Gruszka, A.  and Matthews, G. and Szymura, B. (Eds.) Handbook of Individual Differences in Cognition: Attention, Memory and Executive Control, Springer. }

\seealso{\code{\link{bi.bars}} to show the data by age and gender, \code{\link{irt.fa}} for item factor analysis applying the irt model.} 
\examples{
data(bfi)
describe(bfi)
 keys.list <- list(Agree=c(-1,2:5),Conscientious=c(6:8,-9,-10),Extraversion=c(-11,-12,13:15),Neuroticism=c(16:20),Openness = c(21,-22,23,24,-25))
 keys <- make.keys(28,keys.list,item.labels=colnames(bfi))
 score.items(keys,bfi)
 scores <- score.items(keys[1:27,],bfi[1:27]) #don't score age 
 scores
 
}
\keyword{datasets}
