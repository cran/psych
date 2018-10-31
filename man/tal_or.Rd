\name{Tal_Or}
\alias{Tal_Or}
\alias{Tal.Or}
\alias{pmi}
\alias{tctg}
\docType{data}
\title{
Data set testing causal direction in presumed media influence}
\description{
Nurit Tal-Or, Jonanathan Cohen, Yariv Tasfati, and Albert Gunther (2010) examined the presumed effect of media on other people and change in attitudes.  This data set is from Study 2, and examined the effect of presumed influence of the media upon subsequent actions.  It is used as an example of mediation by Hayes (2013) and for the mediate function. 
}
\usage{data("Tal.Or")}
\format{
  A data frame with 123 observations on the following 6 variables.
  \describe{
    \item{\code{cond}}{Experimental Condition: 0 low media importance, 1 high media importance }
    \item{\code{pmi}}{Presumed media influence (based upon the mean of two items}
    \item{\code{import}}{Importance of the issue }
    \item{\code{reaction}}{Subjects rated agreement about possible reactions to the story (mean of 4 items).}
    \item{\code{gender}}{1 = male, 2 = female }
    \item{\code{age}}{a numeric vector}
  }
}
\details{
Tal-Or et al. (2010) examined the presumed effect of the media in two experimental studies.  These data are from study 2. `... perceptions regarding the influence of a news story about an expected shortage in sugar were manipulated indirectly, by manipulating the perceived exposure to the news story, and behavioral intentions resulting from the story were consequently measured." (p 801). 
}
\source{
The data were downloaded from the webpages of Andrew Hayes (https://www.afhayes.com/public/hayes2018data.zip)  supporting the first and second edition of his book.  The name of the original data set was pmi.  (Gender was recoded to reflect the number of X chromosomes).

The original data are from Nurit Tal-Or, Jonathan Cohen, Yariv Tsfati, and Albert C. Gunther and are used with their kind permission.
}
\references{

Nurit Tal-Or, Jonathan Cohen, Yariv Tsfati and Albert C. Gunther (2010), Testing Causal Direction in the Influence of Presumed Media Influence, Communication Research, 37, 801-824.

Hayes, Andrew F. (2013)  Introduction to mediation, moderation, and conditional process analysis: A regression-based approach.  Guilford Press. 

}
\examples{
data(Tal.Or)
mediate(reaction ~ cond + (pmi), data =Tal.Or,n.iter=50) 

}
\keyword{datasets}
