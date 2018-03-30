\name{blot}
\alias{blot}
\docType{data}
\title{Bond's Logical Operations Test -- BLOT
}
\description{35 items for 150 subjects from Bond's Logical Operations Test.  A good example of Item Response Theory analysis using the Rasch model.  One parameter (Rasch) analysis and two parameter IRT analyses produce somewhat different results. 
}
\usage{data(blot)}
\format{
  A data frame with 150 observations on  35 variables. The BLOT was developed as a paper and pencil test for children to measure Logical Thinking as discussed by Piaget and Inhelder. 
 
}
\details{Bond and Fox apply Rasch modeling to a variety of data sets.  This one, Bond's Logical Operations Test, is used as an example of Rasch modeling for dichotomous items.  In their text (p 56), Bond and Fox report the results using WINSTEPS.  Those results are consistent (up to a scaling parameter) with those found by the rasch function in the ltm package. The WINSTEPS seem to produce difficulty estimates with a mean item difficulty of 0, whereas rasch from ltm has a mean difficulty of -1.52.  In addition, rasch seems to reverse the signs of the difficulty estimates when reporting the coefficients and is effectively reporting "easiness".   

However, when using a two parameter model, one of the items (V12) behaves very differently.

This data set is useful when comparing 1PL, 2PL and 2PN IRT models.  
}
\source{The data are taken (with kind permission from Trevor Bond) from the webpage https://www.winsteps.com/BF3/bondfox3.htm and read using read.fwf.  
}
\references{
T.G. Bond. BLOT:Bond's Logical Operations Test. Townsville, Australia: James Cook Univer- sity. (Original work
  published 1976), 1995.

T. Bond and C. Fox. (2007)  Applying the Rasch model: Fundamental measurement in the human sciences. Lawrence Erlbaum, Mahwah, NJ, US, 2 edition.
}

\seealso{ See also the \code{\link{irt.fa}} and associated plot functions. }
\examples{
data(blot)
#not run
#library(ltm)
#bblot.rasch <- rasch(blot, constraint = cbind(ncol(blot) + 1, 1))  #a 1PL model
#blot.2pl <- ltm(blot~z1)  #a 2PL model
#do the same thing with functions in psych
#blot.fa <- irt.fa(blot)  # a 2PN model
#plot(blot.fa)
}
\keyword{datasets}
