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

\details{Using the data contributed by about 125,000 visitors to the SAPA-project.org website, David Condon has developed a hierarchical framework for assessing personality at two levels. The higher level has the familiar five factors that have been studied extensively in personality research since the 1980s -- Conscientiousness, Agreeableness, Neuroticism, Openness, and Extraversion. The lower level has 27 factors that are considerably more narrow. These were derived based on administrations of about 700 public-domain IPIP items to 3 large samples. Condon describes these scales as being "empirically-derived" because relatively little theory was used to select the number of factors in the hierarchy and the items in the scale for each factor (to be clear, he means relatively little personality theory though he relied on quite a lot of sampling and statistical theory). You can read all about the procedures used to develop this framework in  his book/manual. If you would like to reproduce these analyses, you can download the data files from Dataverse (links are also provided in the manual) and compile this script in R (he used knitR). Instructions are provided in the Preface to the manual.
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
corPlot(sc$corrected,numbers=TRUE,cex=.8,xlas=2)
}
\keyword{datasets}
