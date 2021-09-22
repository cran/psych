\name{GERAS}
\alias{GERAS}
\alias{GERAS.items}
\alias{GERAS.keys}
\alias{GERAS.dictionary}
\alias{GERAS.scales}
\docType{data}
\title{Data from Gruber et al, 2020, Study 2: Gender Related Attributes Survey
}
\description{
Gruber et al. (2020) report on the psychometric properties of a multifaceted Gender Related Attributes Survey. Here are the data from their 3 domains (Personality, Cognition and Activities and Interests from their study 2. Eagly and Revelle (in press) include these data in their review of the power of aggregation. The data  are included here as demonstrations of the \code{\link{cohen.d}} and \code{\link{scatterHist}} functions in the psych package and may be used to show the power of aggregation.
}
\usage{data("GERAS")
#These other objects are included in the file
# data("GERAS.scales")
# data("GERAS.dictionary")
# data("GERAS.items")
# data("GERAS.keys")
    }

\format{
 A data frame with 471 observations on the following 51 variables (selected from the original 93) The code numbers are item numbers from the bigger set.
 \describe{
  \item{\code{V15}}{  reckless}
  \item{\code{V22}}{ willing to take risks}
  \item{\code{V11}}{ courageous}
  \item{\code{V6}}{a  adventurous}
  \item{\code{V19}}{  dominant}
  \item{\code{V14}}{  controlling}
  \item{\code{V20}}{  boastful}
  \item{\code{V21}}{  rational}
  \item{\code{V23}}{ analytical}
  \item{\code{V9}}{   pragmatic}
  \item{\code{V44}}{  to find an address for the first time}
  \item{\code{V45}}{   to find a way again}
  \item{\code{V46}}{  to understand equations}
  \item{\code{V50}}{  to follow directions}
  \item{\code{V51}}{  to understand equations}
  \item{\code{V53}}{    day-to-day calculations}
  \item{\code{V48}}{    to write a computer program}
  \item{\code{V69}}{ paintball}
  \item{\code{V73}}{    driving go-cart}
  \item{\code{V71}}{ drinking beer}
  \item{\code{V68}}{ watching action movies}
  \item{\code{V75}}{ playing cards (poker)}
  \item{\code{V72}}{ watching sports on TV}
  \item{\code{V67}}{  doing certain sports (e.g. soccer, ...)}
  \item{\code{V74}}{    Gym (weightlifting)}
  \item{\code{V27}}{  warm-hearted}
  \item{\code{V28}}{   loving}
  \item{\code{V29}}{   caring}
  \item{\code{V26}}{ compassionate}
  \item{\code{V32}}{  delicate}
  \item{\code{V30}}{   tender}
  \item{\code{V24}}{  familiy-oriented}
  \item{\code{V40}}{  anxious}
  \item{\code{V39}}{  thin-skinned}
  \item{\code{V41}}{  careful}
  \item{\code{V55}}{   to explain foreign words}
  \item{\code{V58}}{to find the right words to express certain content}
  \item{\code{V59}}{ synonyms for a word in order to avoid repetitions}
  \item{\code{V60}}{  to phrase a text}
  \item{\code{V54}}{   remembering events from your own life}
  \item{\code{V63}}{    to notice small changes}
  \item{\code{V57}}{    to remember names and faces}
  \item{\code{V89}}{  shopping}
  \item{\code{V92}}{ gossiping}
  \item{\code{V81}}{   watching a romantic movie}
  \item{\code{V80}}{  talking on the phone with a friend}
  \item{\code{V90}}{  yoga}
  \item{\code{V83}}{    rhythmic gymnastics}
  \item{\code{V84}}{  going for a walk}
  \item{\code{V86}}{  dancing}
  \item{\code{gender}}{gender (M=1 F=2)}
 }
}
\details{These 50 items (+ gender) may be formed into scales using the GERAS.keys
The first 10 items are Male Personality, the next 10  are  Female Personality, then 7 and 7 M and F Cognition, then 8 and 8 M and F Activity items.  The Pers, Cog and Act scales are formed from the M-F scales for the three domains.  M and F are the composites of the Male and then the Female scales.  MF.all is the composite of the  M - F scales.  See the GERAS.keys object for scoring directions.
 

"M.pers" "F.pers" "M.cog"  "F.cog"  "M.act"  "F.act"  "Pers"   "Cog"    "Act"    "M"     
 "F"      "MF.all" "gender"


}
\source{Study 2 data  downloaded from the Open Science Framework https://osf.io/42jhr/ 
Used by kind permission of Freya M. Gruber, Tullia Ortner, and Belinda A. Pletzer.}

\references{

Alice H. Eagly and William Revelle, Understanding the Magnitude of Psychological Differences Between Women and Men Requires Seeing the Forest and the Tree. Perspectives in Psychological Science (in press)

Gruber, Freya M. and Distlberger, Eva and Scherndl, Thomas and Ortner, Tuulia M. and Pletzer, Belinda (2020)   Psychometric properties of the multifaceted Gender-Related Attributes Survey {(GERAS)}  European Journal of Psychological Assessment, 36, (4) 612-623.

}
\examples{
data(GERAS)
GERAS.keys  #show the keys
#show the items from the dictionary
psych::lookupFromKeys(GERAS.keys, GERAS.dictionary[,1,drop=FALSE])


#now, use the GERAS.scales to show a scatterHist  plot showing univariate d and bivariate 
# Mahalanobis D.

psych::scatterHist(F ~ M + gender, data=GERAS.scales, cex.point=.3,smooth=FALSE, 
xlab="Masculine Scale",ylab="Feminine Scale",correl=FALSE, 
d.arrow=TRUE,col=c("red","blue"), bg=c("red","blue"), lwd=4, title="Combined  M and F 
scales",cex.cor=2,cex.arrow=1.25, cex.main=2)







}
\keyword{datasets}
