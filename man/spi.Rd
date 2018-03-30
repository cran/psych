\name{spi}
\alias{spi}
\alias{spi.dictionary}
\alias{spi.keys}
\docType{data}
\title{A sample from the SAPA Personality Inventory including an item dictionary and scoring keys.} 

\description{The SPI (SAPA Personality Inventory) is a set of 135 items primarily selected from International Personality Item Pool (ipip.ori.org). This is an example data set collected using SAPA procedures the sapa-project.org web site.  This data set includes 10 demographic variables as well.  The data set with 4000 observations on 145 variables may be used for examples in scale construction and validation, as well as empirical scale construction to predict multiple criteria.  
}
\usage{data("spi")
data(spi.dictionary)
data(spi.keys)
}
\format{
  A data frame with 4000 observations on the following 145 variables. (The q numbers are the SAPA item numbers).
  \describe{
    \item{\code{age}}{Age in years from 11 -90}
    \item{\code{sex}}{Reported biological sex (coded by X chromosones => 1=Male, 2 = Female)}
    \item{\code{health}}{Self rated health  1-5: poor, fair, good, very good, excellent }
    \item{\code{p1edu}}{Parent 1 education}
    \item{\code{p2edu}}{Parent 2 education}
    \item{\code{education}}{Respondents education: less than 12, HS grad, current univ, 
    some univ, associate degree, college degree, in grad/prof, grad/prof degree }
    \item{\code{wellness}}{Self rated "wellnes" 1-2}
    \item{\code{exer}}{Frequency of exercise: very rarely, < 1/month, < 1/wk,
     1 or 2 times/week, 3-5/wk, > 5 times/week}
    \item{\code{smoke}}{never, not last year, < 1/month, <1/week, 1-3 days/week,
     most days, up to 5 x /day, up to 20 x /day,  > 20x/day}
    \item{\code{ER}}{Emergency room visits none, 1x, 2x, 3 or more times}
    \item{\code{q_253}}{ see the spi.dictionary for these items (q_253}
    \item{\code{q_1328}}{see the dictionary for all items q_1328)}
   }
 
}

\details{Using the data contributed by about 125,000 visitors to the \url{https://SAPA-project.org} website, David Condon has developed a hierarchical framework for assessing personality at two levels. The higher level has the familiar five factors that have been studied extensively in personality research since the 1980s -- Conscientiousness, Agreeableness, Neuroticism, Openness, and Extraversion. The lower level has 27 factors that are considerably more narrow. These were derived based on administrations of about 700 public-domain IPIP items to 3 large samples. Condon describes these scales as being "empirically-derived" because relatively little theory was used to select the number of factors in the hierarchy and the items in the scale for each factor (to be clear, he means relatively little personality theory though he relied on quite a lot of sampling and statistical theory). You can read all about the procedures used to develop this framework in  his book/manual. If you would like to reproduce these analyses, you can download the data files from Dataverse (links are also provided in the manual) and compile this script in R (he used knitR). Instructions are provided in the Preface to the manual.

This small subset of the data is provided for demonstration purposes. 
}
\source{
https://sapa-project.org/research/SPI/SPIdevelopment.pdf.
}
\references{Condon, D. (2017) The SAPA Personality Inventory:An empirically-derived, hierarchically-organized self-report personality assessment model 

}
\examples{
data(spi)
data(sapa.dictionary)
bestScales(spi, criteria="health",dictionary=spi.dictionary)

sc <- scoreItems(spi.keys,spi)
corPlot(sc$corrected,numbers=TRUE,cex=.4,xlas=2,min.length=6)
}
\keyword{datasets}
