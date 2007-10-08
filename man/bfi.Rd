\name{bfi}
\alias{bfi}
\docType{data}
\title{25 Personality items representing 5 factors}
\description{25 personality self report items taken from the International Personality Item Pool (ipip.ori.org) were included as part of the Synthetic Aperture Personality Assessment (SAPA) web based personality assessment project.  The data from 1000 subjects are included here as a demonstration set for scale construction and factor analysis.
}
\usage{data(bfi)}
\format{
  A data frame with 1000 observations on the following 25 variables.
  \describe{
    \item{\code{A1}}{Am indifferent to the feelings of others.}
    \item{\code{A2}}{Inquire about others' well-being.}
    \item{\code{A3}}{Know how to comfort others.}
    \item{\code{A4}}{Love children.}
    \item{\code{A5}}{Make people feel at ease.}
    \item{\code{C1}}{Am exacting in my work.}
    \item{\code{C2}}{Continue until everything is perfect.}
    \item{\code{C3}}{Do things according to a plan.}
    \item{\code{C4}}{Do things in a half-way manner.}
    \item{\code{C5}}{Waste my time.}
    \item{\code{E1}}{Don't talk a lot.}
    \item{\code{E2}}{Find it difficult to approach others.}
    \item{\code{E3}}{Know how to captivate people.}
    \item{\code{E4}}{Make friends easily.}
    \item{\code{E5}}{Take charge.}
    \item{\code{N1}}{Get angry easily.}
    \item{\code{N2}}{Get irritated easily.}
    \item{\code{N3}}{Have frequent mood swings.}
    \item{\code{N4}}{Often feel blue.}
    \item{\code{N5}}{Panic easily.}
    \item{\code{O1}}{Am full of ideas.}
    \item{\code{O2}}{Avoid imposing my will on others.}
    \item{\code{O3}}{Carry the conversation to a higher level.}
    \item{\code{O4}}{Spend time reflecting on things.}
    \item{\code{O5}}{Will not probe deeply into a subject.}
  }
}
\details{The 25 items are organized by five putative factors: Agreeableness, Conscientiousness, Extraversion, Neuroticism, and Opennness.  
}
\source{The items are from the ipip.  The data are from the SAPA project, collected fall,2006.
}
\references{
}
\examples{
data(bfi)
describe(bfi)
keys <- matrix(c(-1,1,1,1,1,rep(0,25),1,1,1,-1,-1,rep(0,25),-1,-1,1,1,1,rep(0,25),1,1,1,1,1,rep(0,25),1,-1,1,1,-1),ncol=5)
rownames(keys) <- colnames(bfi)
colnames(keys) <- c("Agreeable","Conscientious","Extravert","Neurotic","Open")
score.items(keys,bfi,short=TRUE)

}
\keyword{datasets}
