\name{burt}
\alias{burt}
\docType{data}
\title{11 emotional variables from Burt (1915)}

\description{Cyril Burt reported an early factor analysis with a circumplex structure of 11 emotional variables in 1915.  8 of these were subsequently used by Harman in his text on factor analysis.  Unfortunately, it seems as if Burt made a mistake for the matrix is not positive definite.  With one change from .87 to .81 the matrix is positive definite.
}
\usage{data(burt)}
\format{
  A correlation matrix based upon 172 "normal school age children aged 9-12". 
  \describe{
    \item{Sociality}{Sociality}
    \item{Sorrow}{Sorrow}
    \item{Tenderness}{Tenderness}
    \item{Joy}{Joy}
    \item{Wonder}{Wonder}
    \item{Elation}{Elation}
    \item{Disgust}{Disgust}
    \item{Anger}{Anger}
    \item{Sex}{Sex}
    \item{Fear}{Fear}
    \item{Subjection}{Subjection}
    }
}

\details{
The Burt data set is interesting for several reasons.  It seems to be an early example of the organizaton of emotions into an affective circumplex, a subset of it has been used for factor analysis examples (see \code{\link{Harman.Burt}}, and it is an example of how typos affect data.  The original data matrix has one negative eigenvalue. With the replacement of the correlation between Sorrow and Tenderness from .87 to .81, the matrix is positive definite.

Alternatively, using \code{\link{cor.smooth}}, the matrix can be made positive definite as well, although cor.smooth makes more (but smaller) changes.
}
\source{
(retrieved from the web at http://www.biodiversitylibrary.org/item/95822#790)
Following a suggestion by Jan DeLeeuw.

}
\references{
Burt, C.General and Specific Factors underlying the Primary Emotions. Reports of the British Association for the Advancement of Science, 85th meeting, 
held in Manchester,  September 7-11, 1915. 
London, John Murray, 1916, p. 694-696 
(retrieved from the web at http://www.biodiversitylibrary.org/item/95822#790) 

}


\seealso{ \code{\link{Harman.Burt}} in the \code{\link{Harman}} dataset and \code{\link{cor.smooth}}
}

\examples{
data(burt)
eigen(burt)$values  #one is negative!
burt.new <- burt
burt.new[2,3] <- burt.new[3,2] <- .81
eigen(burt.new)$values  #all are positive
bs <- cor.smooth(burt)
round(burt.new - bs,3)

}
\keyword{datasets}
