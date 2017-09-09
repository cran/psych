\name{ability}
\alias{ability}
\docType{data}
\title{16 ability items scored as correct or incorrect.}
\description{
16 multiple choice ability items 1525 subjects taken from the Synthetic Aperture Personality Assessment (SAPA) web based personality assessment project are saved as \code{\link{iqitems}}. Those data are shown as examples of how to score multiple choice tests and analyses of response alternatives.  When scored correct or incorrect, the data are useful for demonstrations of  tetrachoric based factor analysis \code{\link{irt.fa}} and finding tetrachoric correlations. 
}
\usage{data(iqitems)}
\format{
   A data frame with 1525 observations on the following 16 variables. The number following the name is the item number from SAPA.
  \describe{
    \item{\code{reason.4}}{Basic reasoning questions }
    \item{\code{reason.16}}{Basic reasoning question}
    \item{\code{reason.17}}{Basic reasoning question}
    \item{\code{reason.19}}{Basic reasoning question }   
     \item{\code{letter.7}}{In the following alphanumeric series, what letter comes next?}
    \item{\code{letter.33}}{In the following alphanumeric series, what letter comes next?}
    \item{\code{letter.34}}{In the following alphanumeric series, what letter comes next}
    \item{\code{letter.58}}{In the following alphanumeric series, what letter comes next?}
    \item{\code{matrix.45}}{A matrix reasoning task}
    \item{\code{matrix.46}}{A matrix reasoning task}
    \item{\code{matrix.47}}{A matrix reasoning task}
    \item{\code{matrix.55}}{A matrix reasoning task}
    \item{\code{rotate.3}}{Spatial Rotation of type 1.2}
    \item{\code{rotate.4}}{Spatial Rotation of type 1.2}
     \item{\code{rotate.6}}{Spatial Rotation of type 1.1}
    \item{\code{rotate.8}}{Spatial Rotation of type 2.3}
  }
}
\details{16 items were sampled from 80 items given as part of the SAPA (\url{http://sapa-project.org}) project  (Revelle, Wilt and Rosenthal, 2009; Condon and Revelle, 2014) to develop online measures of ability. These 16 items reflect four lower order factors (verbal reasoning, letter series, matrix reasoning, and spatial rotations.  These lower level factors all share a higher level factor ('g').  

This data set may be used to demonstrate item response functions, \code{\link{tetrachoric}} correlations, or \code{\link{irt.fa}} as well as \code{\link{omega}} estimates of of reliability and hierarchical structure.

In addition, the data set is a good example of doing item analysis to examine the empirical response probabilities of each item alternative as a function of the underlying latent trait.  When doing this, it appears that two of the matrix reasoning problems do not have monotonically increasing trace lines for the probability correct.  At moderately high ability (theta = 1) there is a decrease in the probability correct from theta = 0 and theta = 2.
}
\source{ The example data set is taken from the Synthetic Aperture Personality Assessment personality and ability test at \url{http://sapa-project.org}.  The data were collected with David Condon from 8/08/12 to 8/31/12.

Similar data are available from the International Cognitive Ability Resource at \url{http://icar-project.com}.
} 

\references{Revelle, William, Wilt, Joshua,  and Rosenthal, Allen (2010)  Personality and Cognition: The Personality-Cognition Link. In Gruszka, Alexandra  and Matthews, Gerald   and Szymura, Blazej (Eds.) Handbook of Individual Differences in Cognition: Attention, Memory and Executive Control, Springer. 

Condon, David and Revelle, William, (2014) The International Cognitive Ability Resource:  Development and initial validation of a public-domain measure. Intelligence, 43, 52-64.
}
\examples{
data(ability)
#not run
# ability.irt <- irt.fa(ability)
# ability.scores <- score.irt(ability.irt,ability)
}
\keyword{datasets}
